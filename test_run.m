% example run:
ix = 1;
iy = 1;
% dirName = '/export/scratch/albertas/download_dump/S. pyogenes all data';
dirName = '/export/scratch/albertas/download_dump/pyo_test/';
depth=1
windowWidths=[0 250:50:1000];
sF = 0.9:0.025:1.1;
thryFiles = dir('/export/scratch/albertas/download_dump/single/theoryOutput/*.mat');


local_pipeline_mp(ix, iy, dirName, depth, windowWidths, sF, thryFiles)

 c = parcluster;
% c.AdditionalProperties.AccountName = 'snic2021-5-132';
% c.AdditionalProperties.WallTime = '1:00:00';
j =batch(c,@local_pipeline_mp,1,{ix, iy, dirName, depth, 0, sF, thryFiles},'Pool',30);

%%
% example run:
ix=2
iy=1
dirName = '/proj/snic2022-5-384/users/x_albdv/data/discriminative/pyo/';
depth=1
windowWidths=0 ; %[0 250:50:1000];
sF = 0.9:0.025:1.1;
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files April 2024/*.mat');


local_pipeline_mp(ix, iy, dirName, depth, windowWidths, sF, thryFiles)


