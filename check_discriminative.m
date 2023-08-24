file = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local_Pyogenes/sample358/211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=400_table_2023-06-16_01_37_18.txt';
%%
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat');

thryFileIdx = 1;
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);


sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);


import Core.extract_species_name; % find e-coli
[speciesLevel, idc] = extract_species_name(theoryStruct);

%% -
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local_Pyogenes/'))
file = '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt';
% file = '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt';

filePCC = '211213_Sample358-3-st2_110nm_0.169nmPERbp_PCC_table_2023-06-15_23_41_26.txt';
[rezMax,barnames] = Core.load_coefs(filePCC);



import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allNums,refNums,signMatch, fp,positives] = ...
    discrim_true_positives(rezMax, speciesLevel, idc);

positivesPCC = positives==1;
% theoryStruct([cell2mat(refNums(2))]).name

{theoryStruct([cell2mat(refNums(2))]).name}'


filesMP = {'210901_Sample2_110nmPERpx_0.200nmPERbp_MP_w=400_table_2023-06-29_13_18_54.txt'};


filesMP = {
    '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=400_table_2023-06-16_01_37_18.txt',...
    '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt',...
    '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=600_table_2023-06-16_03_43_46.txt'};

successRate = zeros(1,3);
for ix=1:3
    [rezMax,barnamesMP] = Core.load_coefs(filesMP{ix});
%     try
%         barsPassThresh = ismember(barnames,barnamesMP);
%     catch
        barsPassThresh = ismember(cellfun(@(x) strrep(x.name(1:end-4),'-','_'),barGen,'un',false),barnamesMP);
%     end


    [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positivesMP] = ...
        discrim_true_positives(rezMax, speciesLevel, idc);

    positivesMP{ix} = positivesPCC;

    positivesMP{ix}(find(barsPassThresh)) = positives==1;

    successRate(ix) = sum(positivesMP{ix})/length(barnames);
end
% sum(positives==1)/length(rezMax{1})
theoryStruct(refNums{6}).name
{theoryStruct([cell2mat(refNumsMP(25))]).name}'

cellfun(@(x) x.bestBarStretch, rezMax{refNums{6}(1)})


%%
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/'));
file = 'RawKymos_191018_130nm_0.243nmPERbp_MP_w=500_table_2023-06-15_10_24_48.txt';
filePCC = 'RawKymos_191018_130nm_0.243nmPERbp_PCC_table_2023-04-19_16_40_11.txt';


file = 'RawKymos_191018_130nm_0.243nmPERbp_MP_w=400_table_2023-06-15_10_05_21.txt';
filePCC = 'RawKymos_191018_130nm_0.243nmPERbp_PCC_table_2023-04-19_16_40_11.txt'
bar = 'x1_6x_DA65788_191018_Fragment_27_ZVIExport_24_molecule_1_kymogr';


file = 'RawKymographs - 190903_130nm_0.261nmPERbp_MP_w=600_table_2023-06-15_10_37_11.txt';
filePCC = 'RawKymographs - 190903_130nm_0.261nmPERbp_PCC_table_2023-04-19_16_51_32.txt';
bar = 'x1_6x_DA65793_190903_Fragment_2_ZVIExport_02_molecule_1_kymogra'
% file = '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt';
% file = '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt';

% filePCC = '211213_Sample358-3-st2_110nm_0.169nmPERbp_PCC_table_2023-06-15_23_41_26.txt';
[rezMax,barnames] = Core.load_coefs(filePCC);

[barid,cc] = ismember({bar},barnames);


import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allNums,refNums,signMatch, fp,positives] = ...
    discrim_true_positives(rezMax, speciesLevel, idc);

positivesPCC = positives==1;
% theoryStruct([cell2mat(refNums(2))]).name

{theoryStruct([cell2mat(refNums(2))]).name}'


filesMP = {file};


% filesMP = {
%     '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=400_table_2023-06-16_01_37_18.txt',...
%     '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt',...
%     '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=600_table_2023-06-16_03_43_46.txt'};

successRate = zeros(1,3);
for ix=1:length(filesMP)
    [rezMax2,barnamesMP] = Core.load_coefs(filesMP{ix});

    [barid2,cc2] = ismember({bar},barnamesMP);

%     try
%         barsPassThresh = ismember(barnames,barnamesMP);
%     catch
        barsPassThresh = ismember(cellfun(@(x) strrep(x.name(1:end-4),'-','_'),barGen,'un',false),barnamesMP);
%     end


    [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positivesMP] = ...
        discrim_true_positives(rezMax2, speciesLevel, idc);

    positivesMP{ix} = positivesPCC;

    positivesMP{ix}(find(barsPassThresh)) = positives==1;

    successRate(ix) = sum(positivesMP{ix})/length(barnames);
end
% sum(positives==1)/length(rezMax{1})
theoryStruct(refNums{11}).name
{theoryStruct([cell2mat(refNumsMP(2))]).name}'

cellfun(@(x) x.bestBarStretch, rezMax{refNums{6}(1)})


%%

