clear all
close all
% clc
% figure('Menu','none','ToolBar','none');

projectName = 'vesselReg';
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
        %%%% slicer
        src.slicer = 'ml slicer/5.0.3';
        system([src.slicer '; slicer --version > /dev/null'],'-echo');
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

% Skeletonize the vessel segmentation using MATLAB
mriTofSeg = MRIread(tof_vesselSeg);
mriTofSkeleton = mriTofSeg;
mriTofSkeleton.vol = bwskel(mriTofSeg.vol > 0);
tof_skeleton = fullfile(fileparts(tof_vesselSeg), 'tof_skeleton.nii.gz');
% MRIwrite(mriTofSkeleton, tof_skeleton);



% segAndSkeleton = mriTofSeg.vol;
% segAndSkeleton(mriTofSkeleton.vol>0) = 2;
% segAndSkeleton = uint8(segAndSkeleton./max(segAndSkeleton(:)) .* 255);
% orthosliceViewer(segAndSkeleton);


% % Find connected components --- from skeleton
% CC = bwconncomp(mriTofSkeleton.vol, 26);
% [a,b] = sort(cellfun(@length, CC.PixelIdxList),'descend');
% CC.PixelIdxList = CC.PixelIdxList(b);

% tof_skeleton_label = fullfile(fileparts(tof_skeleton), 'tof_skeleton_label.nii.gz');
% mriTofSkeletonLabel = mriTofSkeleton;
% mriTofSkeletonLabel.vol = zeros(size(mriTofSkeleton.vol));
% for i = 1:length(CC.PixelIdxList)
%     mriTofSkeletonLabel.vol(CC.PixelIdxList{i}) = i;
% end
% MRIwrite(mriTofSkeletonLabel, tof_skeleton_label);
% % manually select good components
% okCompIdx_fromSkel = [2 31 21 15 24 19 51 18];



% Find connected components --- from seg
% CC = bwconncomp(mriTofSeg.vol > 0, 26);
% [a,b] = sort(cellfun(@length, CC.PixelIdxList),'descend');
% CC.PixelIdxList = CC.PixelIdxList(b);

tof_seg_label = fullfile(fileparts(tof_vesselSeg), 'tof_seg_label.nii.gz');
mriTofSegLabel = mriTofSeg;
mriTofSegLabel.vol = zeros(size(mriTofSeg.vol));
% for i = 1:length(CC.PixelIdxList)
%     mriTofSegLabel.vol(CC.PixelIdxList{i}) = i;
% end
% MRIwrite(mriTofSegLabel, tof_seg_label);
% manually select good components
okCompIdx_fromSeg = [33 17 19 20 13 52];
% write each good component to a separate nii file
mri = mriTofSeg;
for i = 1:length(okCompIdx_fromSeg)
    mri.vol = zeros(size(mri.vol));
    % mri.vol(CC.PixelIdxList{okCompIdx_fromSeg(i)}) = 1;
    tof_vessel_label{i} = fullfile(fileparts(tof_seg_label), ['tof_vessel_' num2str(okCompIdx_fromSeg(i)) '.nii.gz']);
    % MRIwrite(mri, tof_vessel_label{i});
end

return

forceThis = 1;
for i = 5%1:length(okCompIdx_fromSeg)
    tof_vessel_centerlines{i} = replace(replace(tof_vessel_label{i}, '.nii.gz', '.vtk'), '/seg/', '/centerlines/');
    if ~exist(fileparts(tof_vessel_centerlines{i}),'dir'); mkdir(fileparts(tof_vessel_centerlines{i})); end
    if forceThis || ~exist(tof_vessel_centerlines{i},'file')
        slicer_centerline_from_mask(tof_vessel_label{i}, tof_vessel_centerlines{i});
        % disp(strjoin(cmd,newline));
    end
end
vmtk_viewVolAndSurf(tof, tof_vessel_centerlines{i})





% create surface for each good component from seg (marching cubes -> clean -> smoothing [-> upsample])
doSurfUpSample = true;  % set true to upsample for thin/coarse meshes (butterfly, 1 pass)
for i = 1:length(okCompIdx_fromSeg)
    tof_vessel_surf_files{i} = replace(replace(tof_vessel_label{i}, '.nii.gz', '.vtk'), '/seg/', '/surf/');
    if ~exist(fileparts(tof_vessel_surf_files{i}),'dir'); mkdir(fileparts(tof_vessel_surf_files{i})); end
    if forceThis || ~exist(tof_vessel_surf_files{i},'file');
        tof_vessel_surf_raw    = replace(tof_vessel_surf_files{i}, '.vtk', '_raw.vtk');
        tof_vessel_surf_cleaned = replace(tof_vessel_surf_files{i}, '.vtk', '_cleaned.vtk');
        vmtk_surfFromSeg(tof_vessel_label{i}, tof_vessel_surf_raw);
        vmtk_surfClean(tof_vessel_surf_raw, tof_vessel_surf_cleaned);
        vmtk_surfSmoothing(tof_vessel_surf_cleaned, tof_vessel_surf_files{i});
        if doSurfUpSample
            tof_vessel_surf_up = replace(tof_vessel_surf_files{i}, '.vtk', '_up.vtk');
            vmtk_surfUpSample(tof_vessel_surf_files{i}, tof_vessel_surf_up);
            movefile(tof_vessel_surf_up, tof_vessel_surf_files{i});
        end
        vmtk_viewVolAndSurf(tof_vessel_label{i}, tof_vessel_surf_files{i})
    end
end



for i = 1:length(okCompIdx_fromSeg)
    tof_vessel_centerlines{i} = replace(tof_vessel_surf_files{i}, '/surf/', '/centerlines/');
    if ~exist(fileparts(tof_vessel_centerlines{i}),'dir'); mkdir(fileparts(tof_vessel_centerlines{i})); end
    if forceThis || ~exist(tof_vessel_centerlines{i},'file')
        cmd = vmtk_centerlinesFromSurf(tof_vessel_surf_files{i}, tof_vessel_centerlines{i},[],[],1);
        disp(strjoin(cmd,newline));
    end
end

cmd = vmtk_viewVolAndCenterlines(tof, tof_vessel_centerlines{i}, 1);
disp(strjoin(cmd,newline));











forceThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract vessel centerlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

forceThis = 1;
% Step 1: Create surface mesh from segmentation
tof_surf_file = fullfile(projectScratch,'tof','surf','tof_surface.vtk');
if ~exist(fileparts(tof_surf_file),'dir'); mkdir(fileparts(tof_surf_file)); end
if forceThis || ~exist(tof_surf_file,'file');
    vmtk_surfFromSeg(tof_vesselSeg, tof_surf_file);
end
vmtk_viewVolAndSurf(tof_vesselSeg, tof_surf_file)



% Step 2: Extract centerlines from surface mesh
tof_centerlines_file = fullfile(projectScratch,'tof','centerlines','tof_centerlines.vtk');
if ~exist(fileparts(tof_centerlines_file),'dir'); mkdir(fileparts(tof_centerlines_file)); end
if forceThis || ~exist(tof_centerlines_file,'file');
    vmtk_centerlinesFromSurf(tof_surf_file, tof_centerlines_file, 1);
end
vmtk_viewVolAndCenterlines(tof_vesselSeg, tof_centerlines_file, 1)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Register vessel centerlines vfMRI
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For later.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
