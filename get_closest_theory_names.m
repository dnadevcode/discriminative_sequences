function [] = get_closest_theory_names(dirName)

%% get experiment folders
dirName = '/export/scratch/albertas/data_temp/Alignment/data/ecoli/'; % todo: the code itself takes 5 elements

% thryFiles = dir('/export/scratch/albertas/download_dump/single/theoryOutput/*.mat');
thryFiles = dir('/export/scratch/albertas/data_temp/Alignment/data/theory/*.mat');

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

%%


% here just need the names of theories..
sets.theoryFile{1} = outputLocs{idd}{3};
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);


load(outputLocs{idd}{1}{1})
% batch(c,@run_pipeline_scores,1,{dirName,[twoList(idd,:)], depth, w, sF,timeFramesNr, thryFiles},'Pool',numW);

tic
maxCoef = cell(1,size(matAllCoefs,1));
%% Now find the best coefficient from matAllCoefs (using a cascading or whatever scheme)
for barid =1:size(matAllCoefs,1)
    [singleCoef , singlePos ] =  max(matAllCoefs(barid,:,:),[],2);
    maxCoef{barid} =  squeeze(singleCoef);
end
toc

% unique species names
import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name({theoryStruct.name});

% discriminative locations
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoef}, 0.05);

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);

cellfun(@(x) x(1),refNums)

uniqueSpeciesNames(idSpecies(cellfun(@(x) x(1),refNums(find(is_distinct)))))'

mostCommonSequence{idd} = mode(sort(cell2mat(refNums(find(is_distinct))')));
% uniqueSpeciesNames(idSpecies(refNums{2}))'

end

