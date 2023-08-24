% file = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local_Pyogenes/sample358/211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=400_table_2023-06-16_01_37_18.txt';
%%
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat');

thryFileIdx = 2;
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);


sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);


sets.nmbp = 0.261;
% extract from name
sets.theory.nmbp = sets.nmbp;

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct,sets );



import Core.extract_species_name; % find e-coli
[speciesLevel, idc] = extract_species_name(theoryStruct);

%

%%
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/'));
file = 'RawKymos_191018_130nm_0.243nmPERbp_MP_w=500_table_2023-06-15_10_24_48.txt';
filePCC = 'RawKymos_191018_130nm_0.243nmPERbp_PCC_table_2023-04-19_16_40_11.txt';
bar = 'x1_6x_DA65788_191018_Fragment_25_ZVIExport_23_molecule_1_kymogr';
barcodes = 'RawKymos_191018_130nm_0.243nmPERbpbars.mat';
% 
% 
file = 'RawKymos_191018_130nm_0.243nmPERbp_MP_w=400_table_2023-06-15_10_05_21.txt';
filePCC = 'RawKymos_191018_130nm_0.243nmPERbp_PCC_table_2023-04-19_16_40_11.txt'
bar = 'x1_6x_DA65788_191018_Fragment_27_ZVIExport_24_molecule_1_kymogr';
barcodes = 'RawKymos_191018_130nm_0.243nmPERbpbars.mat';
% 

file = 'RawKymographs - 190903_130nm_0.261nmPERbp_MP_w=600_table_2023-06-15_10_37_11.txt';
filePCC = 'RawKymographs - 190903_130nm_0.261nmPERbp_PCC_table_2023-04-19_16_51_32.txt';
bar = 'x1_6x_DA65793_190903_Fragment_2_ZVIExport_02_molecule_1_kymogra'

% filePCC = '211213_Sample358-3-st2_110nm_0.169nmPERbp_PCC_table_2023-06-15_23_41_26.txt';
[rezMax,barnames] = Core.load_coefs(filePCC);

[barid,cc] = ismember({bar},barnames);
[barsPassThresh,idxkymo] = ismember(cellfun(@(x) matlab.lang.makeValidName(x.name),barGen,'un',false),{bar});

%         N = matlab.lang.makeValidName(barGen{1}.name)

%         N = matlab.lang.makeValidName(barGen{1}.name)

import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allNums,refNums,signMatch, fp,positives] = ...
    discrim_true_positives(rezMax, speciesLevel, idc);

positivesPCC = positives==1;
% theoryStruct([cell2mat(refNums(2))]).name

{theoryStruct([cell2mat(refNums(cc))]).name}'


filesMP = {file};
[rezMax2,barnamesMP] = Core.load_coefs(filesMP{1});
[barid2,cc2] = ismember({bar},barnamesMP);

  [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positivesMP] = ...
        discrim_true_positives(rezMax2, speciesLevel, idc);

{theoryStruct([cell2mat(refNumsMP(cc2))]).name}'


% % successRate = zeros(1,3);
% % for ix=1:length(filesMP)
%     [rezMax2,barnamesMP] = Core.load_coefs(filesMP{ix});
% 
%     [barid2,cc2] = ismember({bar},barnamesMP);
% 
% %     try
% %         barsPassThresh = ismember(barnames,barnamesMP);
% %     catch
%         barsPassThresh = ismember(cellfun(@(x) strrep(x.name(1:end-4),'-','_'),barGen,'un',false),barnamesMP);
% %     end
% 
% 
%     [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positivesMP] = ...
%         discrim_true_positives(rezMax2, speciesLevel, idc);
% 
%     positivesMP{ix} = positivesPCC;
% 
%     positivesMP{ix}(find(barsPassThresh)) = positives==1;
% 
%     successRate(ix) = sum(positivesMP{ix})/length(barnames);
% end
% % sum(positives==1)/length(rezMax{1})
% theoryStruct(refNums{11}).name
% {theoryStruct([cell2mat(refNumsMP(2))]).name}'




% quick_visual_plot(16,9242,barGen,rezMax,bestBarStretch,theoryStruct)
% 
% quick_visual_plot(1,refNumsMP{wIdx}{1}(1),barGen,rezMaxMP,bestBarStretchMP,theoryStruct)


%%
%% Make MP plot simple also work for this!
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = sets.theory.stretchFactors ; % length re-scaling factor
sets.w = 600;
% ixx = 626;

ixx = 2231;%18734;

% quick_visual_plot(1,1,barGen(find(idxkymo)),{{rezMax{ixx}{find(idxkymo)}}},{rezMax{ixx}{find(idxkymo)}.bestBarStretch},theoryStruct(ixx))


bGNew = {};
barGenTest = barGen(find(idxkymo));
bGNew{1} = barGenTest{1};
% bGNew{1}.rawBarcode =  bGNew{1}.rawBarcode(bGNew{1}.rawBitmask);
% bGNew{1}.rawBitmask =  bGNew{1}.rawBitmask(bGNew{1}.rawBitmask);

bGNew{2}.rawBarcode = theoryStruct(ixx).theoryBarcode;
bGNew{2}.rawBitmask = true(1,length(theoryStruct(ixx).theoryBarcode));

%     barcodeGen = bG{idxRun};
    lengths = cellfun(@(x) sum(x.rawBitmask),bGNew);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(bGNew,sF,  sets.w ,timestamp);
    toc
% 
% sets.comparisonMethod = 'mpnan';
% % sets.w = 494;
% 
% [rezMaxMP,bestBarStretchMP,bestLengthMP,rezMaxAllMP] = compare_distance(barGen(find(idxkymo)),theoryStruct(ixx), sets, [] );



import Core.plot_match_simple;
[f] = plot_match_simple(bGNew, oS, 1,2);

%%
% sets.w = 300;
sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
globalov = 0;
% 
% import Core.extract_species_name; % find e-coli
% [speciesLevel, idc] = extract_species_name(theoryStruct);
% if globalov

% compare theory to experiment
import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxT,bestBarStretchT,bestLengthT] = compare_distance(barGen(find(idxkymo)),theoryStruct(ixx), sets, [] );
% )
quick_visual_plot(1,1,barGen(find(idxkymo)),rezMaxT,bestBarStretchT,theoryStruct(ixx))
