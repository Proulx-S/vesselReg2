clear all
close all
% clc
% figure('Menu','none','ToolBar','none');

projectName = 'vesselReg2';
%%%%%%%%%%%%%%%%%%%%%
%% Set up environment
%%%%%%%%%%%%%%%%%%%%%

% Detect computing environment
os   = char(java.lang.System.getProperty('os.name'));
host = char(java.net.InetAddress.getLocalHost.getHostName);
user = char(java.lang.System.getProperty('user.name'));

% Setup folders
if strcmp(os,'Linux') && strcmp(host,'takoyaki') && strcmp(user,'sebp')
    envId      = 1;
    storageDrive = '/local/users/Proulx-S/';
    scratchDrive = '/scratch/users/Proulx-S/';
    projectCode    = fullfile(scratchDrive, projectName);        if ~exist(projectCode,'dir');    mkdir(projectCode);    end
    projectStorage = fullfile(storageDrive, projectName);        if ~exist(projectStorage,'dir'); mkdir(projectStorage); end
    projectScratch = fullfile(scratchDrive, projectName, 'tmp'); if ~exist(projectScratch,'dir'); mkdir(projectScratch); end
    toolDir        = '/scratch/users/Proulx-S/tools';            if ~exist(toolDir,'dir');        mkdir(toolDir);        end
else
    % envId = 2;
    % storageDrive = '/Users/sebastienproulx/bassWrap-reg';
    % scratchDrive = '/Users/sebastienproulx/bassWrap-reg';
    % projectCode    = fullfile(scratchDrive, projectName);        if ~exist(projectCode,'dir');    mkdir(projectCode);    end
    % projectStorage = fullfile(storageDrive, projectName);        if ~exist(projectStorage,'dir'); mkdir(projectStorage); end
    % projectScratch = fullfile(scratchDrive, projectName, 'tmp'); if ~exist(projectScratch,'dir'); mkdir(projectScratch); end
    % toolDir        = '/Users/sebastienproulx/tools';             if ~exist(toolDir,'dir');        mkdir(toolDir);        end
end

% Load dependencies and set paths
%%% matlab util (contains my matlab git wrapper)
tool = 'util'; toolURL = 'https://github.com/Proulx-S/util.git';
if ~exist(fullfile(toolDir, tool), 'dir'); system(['git clone ' toolURL ' ' fullfile(toolDir, tool)]); end; addpath(genpath(fullfile(toolDir,tool)))
%%% matlab others
tool = 'freesurfer'; subTool = 'matlab'; repoURL = 'https://github.com/freesurfer/freesurfer.git';
gitClone(repoURL, fullfile(toolDir, tool), subTool);
%%% bassWrap-reg tools
tool = 'bassWrap-reg';
if exist(fullfile(toolDir, tool), 'dir'); addpath(fullfile(toolDir, tool)); end

%%% neurodesk
switch envId
    case 1
        global src    
        setenv('SINGULARITY_BINDPATH',strjoin({projectCode projectStorage projectScratch toolDir},','));
        %%%% vesselboost
        src.vesselboost = 'ml vesselboost/1.0.0';
        system([src.vesselboost '; prediction.py --help > /dev/null'],'-echo');
        vesselBoostModel = fullfile(projectScratch,'manual_0429');
        system([src.vesselboost '; osf -p abk4p fetch /pretrained_models/manual_0429 ' vesselBoostModel],'-echo');
        %%%% ants
        src.ants = 'ml ants/2.5.3';
        system([src.ants '; antsRegistration --version > /dev/null'],'-echo');
        %%%% freesurfer
        src.fs   = 'ml freesurfer/8.0.0';
        system([src.fs   '; mri_convert > /dev/null'],'-echo');
        %%%% vmtk
        src.vmtk = 'ml vmtk/1.5.0';
        system([src.vmtk '; vmtkcenterlines --help > /dev/null'],'-echo');
        %%%% nipype
        src.nipype = 'ml nipype/1.8.3';
        system([src.nipype '; python -c "import nibabel, skimage.morphology; print(''nipype OK'')"'],'-echo');
    case 2
        warning('neurodesk not implemented for this environment');
    otherwise
        dbstack; error('not implemented')
        % neurodeskModule = {
        % ":/neurodesktop-storage/containers/freesurfer_8.0.0_20250210"
        % ":/neurodesktop-storage/containers/afni_24.3.00_20241003"};
        % for i = 1:length(neurodeskModule)
        %     if contains(getenv("PATH"),neurodeskModule{i}); continue; end
        %     setenv("PATH",getenv("PATH") + neurodeskModule{i});
        % end
end
%% %%%%%%%%%%%%%%%%%%
disp(projectCode)
disp(projectStorage)
disp(projectScratch)


if 0
    %%%%%%%%%%%%%%%%%%
    %% Copy data files
    %%%%%%%%%%%%%%%%%
    % get data pointer
    tmp = fullfile(storageDrive,'vsmCenSur','dataPointers.mat');
    pointerFile = fullfile(projectScratch,'dataPointers.mat');
    copyfile(tmp,pointerFile);

    % get nifis
    load(pointerFile);
    s = 1; % perferct all over
    s = 2; % somewhate rigid movement mostly in the first and 2nd run -> mostly correctable with matlab registration over vesselRegion, but frames failed (will probably require temporal smoothing)
    s = 3; % ok
    s = 4; % minimal possibly non-rigid movement
    s = 5; % some non-rigid movement particularly in one vessel on the left on the last run
    s = 6; % ok
    s = 7; % some movements, not clear if rigid

    S=2;
    %tof
    in = fullfile(roi{S}.(acq).(task).rCond.volAnat.tof.folder,roi{S}.(acq).(task).rCond.volAnat.tof.name);
    out = fullfile(projectCode,'data','tof'); if ~exist(out,'dir'); mkdir(out); end
    tof = fullfile(out,'tof.nii.gz');
    copyfile(in,tof);

    %vfMRI
    in = roi{S}.(acq).(task).rCond.fOrigList{1};
    out = fullfile(projectCode,'data','vfMRI'); if ~exist(out,'dir'); mkdir(out); end
    vfMRI = fullfile(out,'vfMRI.nii.gz');
    copyfile(in,vfMRI);
    %% %%%%%%%%%%%%%%%
else
    %%%%%%%%%%%%%%%%%
    %% Get data files
    tof   = fullfile(projectCode,'data','tof'  ,'tof.nii.gz'  );
    vfMRI = fullfile(projectCode,'data','vfMRI','vfMRI.nii.gz');
    %% %%%%%%%%%%%%%%
end


forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vesselboost vessel segmentation of tof volume data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tof_vesselSeg = fullfile(projectScratch,'tof','seg'); if ~exist(tof_vesselSeg,'dir'); mkdir(tof_vesselSeg); end
tof_vesselSeg = fullfile(tof_vesselSeg,'tof.nii.gz');
if forceThis || ~exist(tof_vesselSeg,'file');
    vesselboost_prediction(fileparts(tof),fileparts(tof_vesselSeg),vesselBoostModel,4);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Skeletonize using scikit-image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use scikit-image's skeletonize (Lee method for 3D)
tof_skeleton_skimage = fullfile(fileparts(tof_vesselSeg), 'tof_skeleton_skimage.nii.gz');
skeletonize_nifti(tof_vesselSeg, tof_skeleton_skimage, 'lee', forceThis);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


forceThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Label connected components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tof_skeleton_skimage_label = fullfile(fileparts(tof_skeleton_skimage), 'tof_skeleton_label.nii.gz');
if forceThis || ~exist(tof_skeleton_skimage_label,'file')
    label_connected_components_nifti(tof_skeleton_skimage, tof_skeleton_skimage_label, 3);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%










mri_tof_skeleton_skimage_label = MRIread(tof_skeleton_skimage_label);

% Get component information for compatibility with existing code
mriTofSkeleton = MRIread(tof_skeleton_skimage);
% Note: Components in labeled image are already sorted (1 = largest, 2 = second largest, etc.)
% Get number of components from labeled image
num_components = max(mri_tof_skeleton_skimage_label.vol(:));
% Create CC structure for compatibility (though we'll use labeled image directly)
CC = bwconncomp(mri_tof_skeleton_skimage_label.vol > 0, 26);
% Rebuild PixelIdxList from labeled image to match sorted order
CC.PixelIdxList = cell(1, num_components);
for i = 1:num_components
    CC.PixelIdxList{i} = find(mri_tof_skeleton_skimage_label.vol == i);
end


% manually select good components
okCompIdx_fromSkel = [21 31 19 15 51 9 18];
% okCompIdx_fromSkel = [2 31 21 15 24 19 51 18];
% write each good component to a separate nii file
mri = mriTofSkeleton;
for i = 1:length(okCompIdx_fromSkel)
    mri.vol = zeros(size(mri.vol));
    mri.vol(CC.PixelIdxList{okCompIdx_fromSkel(i)}) = 1;
    tof_skeleton_labels{i} = fullfile(fileparts(tof_skeleton_skimage_label), ['tof_skeleton_' num2str(okCompIdx_fromSkel(i)) '.nii.gz']);
    MRIwrite(mri, tof_skeleton_labels{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract centerlines from skeleton for each vessel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create centerline VTK files for each identified vessel skeleton
tof_skeleton_centerlines = cell(length(okCompIdx_fromSkel), 1);
for i = 1:length(okCompIdx_fromSkel)
    % Extract skeleton voxel indices for this vessel
    voxel_indices = CC.PixelIdxList{okCompIdx_fromSkel(i)};

    tof_skeleton_labels{i}
    
    % Convert to (x,y,z) coordinates (1-based indexing from ind2sub)
    [x, y, z] = ind2sub(size(mriTofSkeleton.vol), voxel_indices);
    voxel_coords = [x, y, z];
    
    % Convert to RAS coordinates
    ras_coords = voxel_to_ras_coords(voxel_coords, mriTofSkeleton);
    
    % Write VTK centerline file
    output_vtk = fullfile(fileparts(tof_skeleton_skimage_label), ['tof_skeleton_centerline_' num2str(okCompIdx_fromSkel(i)) '.vtk']);
    tof_skeleton_centerlines{i} = ras_coords_to_vtk(ras_coords, output_vtk);
    fprintf('Created centerline: %s (%d points)\n', tof_skeleton_centerlines{i}, size(ras_coords, 1));
end





% Find connected components --- from seg
CC = bwconncomp(mriTofSeg.vol > 0, 26);
[a,b] = sort(cellfun(@length, CC.PixelIdxList),'descend');
CC.PixelIdxList = CC.PixelIdxList(b);

tof_seg_label = fullfile(fileparts(tof_vesselSeg), 'tof_seg_label.nii.gz');
mriTofSegLabel = mriTofSeg;
mriTofSegLabel.vol = zeros(size(mriTofSeg.vol));
for i = 1:length(CC.PixelIdxList)
    mriTofSegLabel.vol(CC.PixelIdxList{i}) = i;
end
MRIwrite(mriTofSegLabel, tof_seg_label);
% manually select good components
okCompIdx_fromSeg = [33 17 19 20 13 52];
% write each good component to a separate nii file
mri = mriTofSeg;
for i = 1:length(okCompIdx_fromSeg)
    mri.vol = zeros(size(mri.vol));
    mri.vol(CC.PixelIdxList{okCompIdx_fromSeg(i)}) = 1;
    tof_vessel_labels{i} = fullfile(fileparts(tof_seg_label), ['tof_vessel_' num2str(okCompIdx_fromSeg(i)) '.nii.gz']);
    MRIwrite(mri, tof_vessel_labels{i});
end

return


vmtk_viewVolAndSurf(tof, tof_vessel_centerlines{i})


