%% On Kurt:
% dirName with kymo's
c = parcluster;
dirName = '/export/scratch/albertas/data_temp/Alignment/data/ehec/EHEC data for local alignment/';
dirName = '/export/scratch/albertas/download_dump/S. pyogenes all data/onlyKymo/';
% folder with theory barcodes 
thryFiles = dir('/export/scratch/albertas/download_dump/single/theoryOutput/*.mat');
numW = 29;
saveDir = '/export/scratch/albertas/output_dump/outputPCC2/';
mkdir(saveDir);


%% On tetralith
c = parcluster;
c.AdditionalProperties.AccountName = 'naiss2024-22-957';
c.AdditionalProperties.WallTime = '6:00:00';
numW = 120;


addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/hca'));
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/discriminative_sequences'));
dirName = '/proj/snic2022-5-384/users/x_albdv/data/discriminative/pyo/';
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files April 2024/*.mat');

%% Parameters

% window widths to run
windowWidths=[0];
% length re-scaling factors to run
sF = 0.9:0.025:1.1;
timeFramesNr = 0;
depth = 1; % maybe don't use

%%
import Helper.get_all_folders;
[barN, twoList,foldname] = get_all_folders(dirName);

idd = 1;
% thryId = 1;
% saveDir = 'output';
% outputLocs = [];
name = [saveDir,foldname{twoList(idd,1)}{twoList(idd,2)}];
% j = batch(c,@run_pipeline_scores,2,...
%     {dirName,[twoList(idd,:)], depth, windowWidths, sF,timeFramesNr, thryFiles,name,thryId},...
%     'Pool',numW);
for idd=1:size(twoList,1)
    name = [saveDir,num2str(idd),'_',foldname{twoList(idd,1)}{twoList(idd,2)}];
    [t,outputLocs{idd}] = run_pipeline_scores(dirName,[twoList(idd,:)], depth, windowWidths, sF, timeFramesNr, thryFiles, name,thryId,theoryStruct,bG{idd});
end
%%
j = cell(1,size(twoList,1));
for idd=1:size(twoList,1)
    name = [saveDir,foldname{twoList(idd,1)}{twoList(idd,2)}];

    j{idd} = batch(c,@run_pipeline_scores,2,...
        {dirName,[twoList(idd,:)], depth, windowWidths, sF,timeFramesNr, thryFiles,name,thryId,theoryStruct},...
    'Pool',numW);
end
%%




