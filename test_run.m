% example run:
ix=2
iy=1
dirName = '/export/scratch/albertas/download_dump/S. pyogenes all data';
depth=1
windowWidths=[0 250:50:1000];
sF = 0.9:0.025:1.1;
thryFiles = dir('/export/scratch/albertas/download_dump/single/theoryOutput/*.mat');


local_pipeline_mp(ix, iy, dirName, depth, windowWidths, sF, thryFiles)

