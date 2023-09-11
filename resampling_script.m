% resampling script

pathMain = 'C:\Users\Lenovo\git\dnadevcode\'; % path to reps.

foldToRun = '';

addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))

addpath('/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/')
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/'))


dirName = 'C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment';
import Core.load_data_fun;
[kymoStructs,barGen,sets] = load_data_fun(dirName,1);


nmbp = sets.nmbp;

%% theory loading

import Core.load_theory_structure;
thryFileIdx = 2; % todo: pass directly the theory file here
[theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);

%%

barGenRun = barGen(1);%(1);
w = [];
[rezMax, bestBarStretch, bestLength, rezOut] = local_alignment_assembly(theoryStruct, barGenRun,w);


import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStruct);
% 
idx = 1;
import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);

{theoryStruct([refNums{1}]).name}'


%%

foldCalc = 'C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment\EC6_Ecoli_DA65808';
foldCalc = 'C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment\EC15_Ecoli_DA65783';

files = dir(fullfile(foldCalc,'*.dat'));

% mpid = zeros(1,length(files));
for i=1:length(files)
    f1= strsplit(files(i).name,'w=');
    if length(f1) == 1
        mppc = i;
    else
        f2 = strsplit(f1{2},'_');
        mpval(i) = str2num(f2{1});
    end
end


for i=1:length(files)
    [rM{i},bnames{i}] = Core.load_coefs(fullfile(files(i).folder,files(i).name));
end

% extract how many scores for each barcode..    
import Core.discrim_true_positives;
pos = zeros(1,length(files));
for i=1:length(files)
    [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positivesMP] = ...
        discrim_true_positives(rM{i}, speciesLevel, idc);
    pos(i) = length(cell2mat(refNumsMP));
end


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
addpath(genpath('C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment\'));
file = 'RawKymos_191009_130nm_0.265nmPERbp_MP_w=400_table_2023-06-15_10_32_59.txt';
filePCC = 'RawKymos_191009_130nm_0.265nmPERbp_PCC_table_2023-04-19_11_21_14.txt';

% 
% file = 'RawKymos_191018_130nm_0.243nmPERbp_MP_w=400_table_2023-06-15_10_05_21.txt';
% filePCC = 'RawKymos_191018_130nm_0.243nmPERbp_PCC_table_2023-04-19_16_40_11.txt'
% bar = 'x1_6x_DA65788_191018_Fragment_27_ZVIExport_24_molecule_1_kymogr';
% 
% 
% file = 'RawKymographs - 190903_130nm_0.261nmPERbp_MP_w=600_table_2023-06-15_10_37_11.txt';
% filePCC = 'RawKymographs - 190903_130nm_0.261nmPERbp_PCC_table_2023-04-19_16_51_32.txt';
% bar = 'x1_6x_DA65793_190903_Fragment_2_ZVIExport_02_molecule_1_kymogra'
% file = '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt';
% file = '211213_Sample358-3-st2_110nm_0.169nmPERbp_MP_w=500_table_2023-06-16_02_53_59.txt';

% filePCC = '211213_Sample358-3-st2_110nm_0.169nmPERbp_PCC_table_2023-06-15_23_41_26.txt';
[rezMaxPrec,barnames] = Core.load_coefs(filePCC);

[barid,cc] = ismember({bar},barnames);


import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allNums,refNums,signMatch, fp,positives] = ...
    discrim_true_positives(rezMaxPrec, speciesLevel, idc);

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



