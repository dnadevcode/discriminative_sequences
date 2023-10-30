
pathMain = '/proj/snic2022-5-384/users/x_albdv/reps/';% reps.

foldToRun = '';

addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))


% curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local_Pyogenes/sample365/20220121_Sample 365-st1_110nm_0.225nmPERbp/';

curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/test/sample365/20220121_Sample365-st1_110nm_0.225nmPERbp/';
import Core.load_local_alignment_results_from_files;
[rM, bnames, mpval,thryNames] = load_local_alignment_results_from_files(curDirKymos ); 





import Core.load_theory_structure;
% thryFileIdx = 1; % todo: pass directly the theory file here
% [theoryStruct,sets] = load_theory_structure(0.225,thryFileIdx);


ix = 1;
import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name([],{'Streptococcus pyogenes'},thryNames{ix}');
%

import Core.disc_locs;
[refNums, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{ix});


    % calcs ref and allnums again..
import Core.discrim_true_positives;
 [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positives] = ...
    discrim_true_positives(rM{ix}, speciesLevel, idc);

thryNames{ix}([cell2mat(refNumsMP(2))])



[length(refNums{1}) positives(1)]

%%
thryNames

import Core.Default.read_default_sets;
hcaSets = read_default_sets('shrinksortersets.txt');

import Core.shrink_finder_fun;
[kymoStructsUpdated,kymoKeep] = shrink_finder_fun( hcaSets, kymoStructs, 0)
