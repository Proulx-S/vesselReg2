cd /scratch/users/Proulx-S/vesselReg/
i=20;

vmtkcenterlines \
-ifile tmp/tof/surf/tof_vessel_${i}.vtk \
-ofile tmp/tof/centerlines/tof_seg_label_${i}.vtk



vmtkimagereader \
-ifile data/tof/tof.nii.gz \
--pipe vmtksurfacereader \
-ifile tmp/tof/centerlines/tof_seg_label_${i}.vtk \
--pipe vmtkrenderer \
--pipe vmtkimageviewer   -i @vmtkimagereader.o   -display 0 \
--pipe vmtksurfaceviewer -i @vmtksurfacereader.o -display 1




cp tmp/tof/centerlines/tof_seg_label_${i}.vtk ../vesselReg2/data
cp tmp/tof/centerlines/tof_vessel_${i}.vtk ../vesselReg2/data