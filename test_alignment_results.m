% test alignment results
%
%   Load any alignment file.


pathMain = '/proj/snic2022-5-384/users/x_albdv/reps/';% reps.

foldToRun = '';

addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))


% curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local_Pyogenes/sample365/20220121_Sample 365-st1_110nm_0.225nmPERbp/';

curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/test/sample365/20220121_Sample365-st1_110nm_0.225nmPERbp/';

import Core.load_local_alignment_results_from_files;
[rM, bnames, mpval,thryNames,files] = load_local_alignment_results_from_files(curDirKymos ); 


import Core.load_theory_structure;
thryFileIdx = 1; % todo: only pass single
[theoryStruct,sets] = load_theory_structure(0.17,thryFileIdx);


% calc stouffer scores
[rM] = stouffer_zscore(rM,theoryStruct,mpval,0.07);
f= figure,histogram(cellfun(@(x) x{1}.stoufferScore,rM{2}))

% Alternative: use MAT Files
%% need barGen
scores = cell(1,length(barGen));

for ii=1:length(barGen)
    [scores{ii},pccScore] = local_bootstrap_run( barGen(ii),rM,bnames,theoryStruct ,mpval,speciesLevel,idc);

end
% [scores,pccScore] = local_bootstrap_run( barGen(ii),rM,bnames,theoryStruct ,mpval,speciesLevel,idc);


for ii=1:length(barGen)
    scores{ii}
    
    % output:
    import Core.export_coefs_resampling;
    T = export_coefs_resampling(scores{ii}, barGen(ii), mpval, [curDirKymos, '_resampling_table'],timestamp);

    save( [curDirKymos, '_resampling_table.mat'],'T');
end
%%


% import Core.load_theory_structure;
% thryFileIdx = 1; % todo: pass directly the theory file here
% [theoryStruct,sets] = load_theory_structure(0.225,thryFileIdx);


ix = 1;
import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name([],{'Streptococcus pyogenes'},thryNames{ix}');
% [speciesLevel,idc] = extract_species_name([],{'Escherichia coli'},thryNames{ix}');

%

% import Core.Discriminative.extract_species_name;
% [uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name(thryNames);


import Core.disc_locs;
[refNums, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{ix});
% thryNames{ix}([cell2mat(refNums(11))])


    % calcs ref and allnums again..
import Core.discrim_true_positives;
 [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positives] = ...
    discrim_true_positives(rM{ix}, speciesLevel, idc);

thryNames{ix}([cell2mat(refNumsMP(2))])



[length(refNums{1}) positives(1)]

%% Alternative: new, creates an output table

import Core.Default.read_default_sets;
hcaSets = read_default_sets('hcaalignmentsets.txt');
hcaSets.default.kymofolder{1} = 'kymo';
hcaSets.default.timestamp ='test';
% 
import Core.identify_discriminative;
[is_distinct,numMatchingSpecies,uniqueMatchSequences,refNums,refNumBad] =...
    identify_discriminative(hcaSets.default,barGen,rezMaxMP, thryNames{1});
    

thryNames

%% Try to do resampling plot as well

% resampling_script_fun('/export/scratch/albertas/data_temp/bargrouping/local_Test/',1)




test_new_scores

%% plot

%%
import CBT.Hca.UI.Helper.plot_any_bar;
ix = 1;
thrIx = 1;
plot_any_bar(ix,barGen,rezMaxMP,theoryStruct,refNums{ix}(thrIx));
%%


% shrink?
% import Core.shrink_finder_fun;
% [kymoStructsUpdated,kymoKeep] = shrink_finder_fun( hcaSets, kymoStructs, 0)



%% Check alignment results from stouffer/resampled stuff

% isdiscriminative without bargen
import Core.Default.read_default_sets;
hcaSets = read_default_sets('hcaalignmentsets.txt');
hcaSets.default.kymofolder{1} = 'kymo';
hcaSets.default.timestamp ='test';
% 
is_distinct = cell(1,length(rM));
numMatchingSpecies = cell(1,length(rM));
uniqueMatchSequences = cell(1,length(rM));

allbars =  bnames{find(mpval==0)};

for i=1:length(rM)

    for j=1:length(bnames{i})

        

    end

%     passingThreshBars = bnames{i};
%     if length(passingThreshBars)~= length(rM{j}{1})
%             passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen)*bars.sets.theory.stretchFactors(end) >= mpval(j));
%     end
%     if length(passingThreshBars)~= length(rM{j}{1})
%         passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen) >=  mpval(j));
%     end


    import Core.identify_discriminative;
    [is_distinct{i},numMatchingSpecies{i},uniqueMatchSequences{i},refNums,refNumBad] =...
        identify_discriminative(hcaSets.default,[],rM{i}, thryNames{i});
end

%% connect to correct bar
    
