addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/lldev'));addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/hca'));addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/bargroupingprototype'));addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/discriminative_sequences'));
load('ecoli_2_all.mat');

load('ecoli_2_fig_individual.mat'); 

%
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat');

thryFileIdx = 1;
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);


sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);


sets.nmbp = 0.25;
% extract from name
sets.theory.nmbp = sets.nmbp;

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct,sets );



import Core.extract_species_name; % find e-coli
[speciesLevel, idc] = extract_species_name(theoryStruct);


%


% theoryStruct([refNumsMP{5}]).name;
windowWidths = 400:100:600;
sets.comparisonMethod = 'mpnan';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
sets.theory.stretchFactors = 0.9:0.025:1.1; %as per 


barGen = bG{1}(1:10);

import CBT.Hca.Core.Comparison.compare_distance;

for wIdx = 1:length(windowWidths)
    sets.w = windowWidths(wIdx);
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),barGen) >= sets.w);

    % assign standard scores
%     rezMaxMP = rezMax;
%     bestBarStretchMP = bestBarStretch;
%     bestLengthMP = bestLength;

    [rezMaxMP,bestBarStretchMP,bestLengthMP,rezMaxAllMP] = compare_distance(barGen(passingThreshBars),theoryStruct, sets, [] );

    import Core.discrim_true_positives;
    [truePositivesMP{wIdx},discSpeciesMP{wIdx},discAllMP{wIdx},allNumsMP{wIdx},refNumsMP{wIdx},signMatchMP{wIdx},fp,positives] =...
        discrim_true_positives(rezMaxMP,speciesLevel,idc);

end

% Max position to more full ovelap

% use plot simple!


quick_visual_plot(16,9242,barGen,rezMax,bestBarStretch,theoryStruct)

