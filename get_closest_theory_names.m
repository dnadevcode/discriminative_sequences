function [] = get_closest_theory_names(dirName)

%% get experiment folders
dirName = '/export/scratch/albertas/data_temp/Alignment/data/ecoli_small/'; % todo: the code itself takes 5 elements

thryFiles = dir('/export/scratch/albertas/download_dump/single/theoryOutput/*.mat');

import Helper.get_all_folders;
[barN, twoList] = get_all_folders(dirName);

sF = 1;
w=0;
numW = 30;
timeFramesNr = 0;
depth = 1;

idd = 1;
outputLocs = [];
[t,outputLocs{idd}] = run_pipeline_scores(dirName,[twoList(idd,:)], depth, w, sF,timeFramesNr, thryFiles);


% batch(c,@run_pipeline_scores,1,{dirName,[twoList(idd,:)], depth, w, sF,timeFramesNr, thryFiles},'Pool',numW);


end

