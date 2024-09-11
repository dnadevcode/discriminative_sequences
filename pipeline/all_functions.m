
%% load_kymo_data_from_fold
%   Args:
%         dirName
%         refsNames
%         allTheoryFold
%   uses default settings such as sF factor (0.8:0.025:1.2), minimum 20
%   timeframes
%   Returns:
%       kymoStructs
%       barN - how many barcodes in each subfolder
%       twoList - list of folders and subfolders to run
%       bG - barcodeGen in each subfolder
%       expPar - experimental parameters extracted from each folder name
%       fastaFileF - name of the theory closest matching the experiment
dirName = '/export/scratch/albertas/data_temp/Alignment/data/ehec/EHEC data for local alignment/';
dirName = '/export/scratch/albertas/download_dump/S. pyogenes all data/onlyKymo/';

refsNames = '/home/avesta/albertas/reps/discriminative_sequences/constants/pyo_theories.txt';
allTheoryFold = '/export/scratch/albertas/download_dump/single/*.fasta';

sF = 0.8:0.025:1.2;
[kymoStructs,barN,twoList,bG,expPar,fastaFileF] = load_kymo_data_from_fold(dirName, refsNames,allTheoryFold,sF,0);

%% Load all names and unique seq 
matFile = '/export/scratch/albertas/download_dump/single/theoryOutput/theoryGen_0.34_110_300_0_2024-04-24_19_10_40_session.mat';
[uniqueSpeciesNames,idSpecies,thryNames,theoryStruct] = thryNames_from_mat(matFile);

%% Find all sequences belonging to specific species
fastaFold = '/export/scratch/albertas/download_dump/single/*.fasta';
    seqNameToTest = 'Escherichia coli';
% seqNameToTest = 'Streptococcus pyogenes';
[thryId] = find_specific_seq(fastaFold, seqNameToTest, matFile);

%% Run alignment
run_alignment

%% Merge outputs from the sample samples
newDir = '/export/scratch/albertas/output_dump/outputPCCGrouped/';
mkdir(newDir);
merge_alignment_outputs(saveDir,newDir, twoList)


%% Get closest matching theory
addpath(genpath('/proj/dnadevdata/reps/discriminative_sequences'))
cdiff = 0.05; % c diff
sflevel = 1; % 1 :-10.. 2: -7.5.. 3: -5 .. ...
[mostCommonSequence,mostCommonRep,namesDir,thryDiscName] = get_best_theory(newDir, [], idSpecies(thryId),barN, twoList,foldname,thryNames(thryId),'hca', sflevel,cdiff);



%% [t,outputLocs] = run_pipeline_scores(dirName, selIdxs, depth, windowWidths, sF, timeFramesNr, thryFiles,savedir)
% main function to run pipelinescores


%% run_mp_based_local_alignment
%   Script to run matrix profile based local alignment


%% Pre-generate theory
CBT.SimpleTwoState.pregenerate_ligand_index(sequence_list,'/proj/dnadevdata/users/x_albdv/data/all/single2/',4,1,110/0.34)
