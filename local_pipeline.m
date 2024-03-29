function [] = local_pipeline(ix,dirName,depth,windowWidths)


% EC17 = DA65780
% EC15 = DA65783
% EC6 = DA65808
% EC24 = DA65788
% EC25 = DA65793

if nargin < 4
    windowWidths = 400:100:600;
end

if nargin < 3
    depth = 1;
end
% folds
% directory name
thryFiles = [dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat'),dir('/export/scratch/albertas/data_temp/bargrouping/New Ref-Theory files May 2022/*.mat')];

if nargin < 2
    dirName = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/';
end
addpath(genpath(dirName));
[~,~] =mkdir('output');
if nargin < 1
    ix = 1;
end

dirStruct = dir(dirName);
dirStruct(~[dirStruct.isdir]) = [];  %remove non-directories
dirStruct(ismember( {dirStruct.name}, {'.', '..'})) = [];  %remove . and ..


if depth == 1
    subDir = dir(fullfile(dirStruct(ix).folder,dirStruct(ix).name));
    subDir(ismember( {subDir.name}, {'.', '..'})) = [];  %remove . and ..
    subDir(find(~cell2mat({subDir.isdir}))) = [];
else
    subDir = dirStruct(ix);
end

%
iy = 1; % most likely single run
sets.dirName = fullfile(subDir(iy).folder,subDir(iy).name);

spltName = strsplit(sets.dirName ,'_');
spltName2 = strsplit(spltName{end},'nm');
spltName3 = strsplit(spltName{end-1},'nm');
sets.nmbp = str2double(spltName2{1});
nmpx = str2double(spltName3{1});


% load theory
thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,spltName{end-1}(1:5))),1:length(thryFiles)));
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);

%%

files = dir(fullfile(sets.dirName,'*.tif'));

sets.kymosets.filenames = arrayfun(@(x) files(x).name,1:length(files),'un',false);
sets.kymosets.kymofilefold = arrayfun(@(x) files(x).folder,1:length(files),'un',false);

% simdata = 0;
sets.output.matDirpath = 'output';
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 6;
sets.minLen = 150;
sets.genConsensus  = 0;

%  following "Strain-level bacterial typing directly from patient
% samples using optical DNA mapping"
sets.timeFramesNr = 20;
sets.theory.stretchFactors = 0.9:0.025:1.1; %as per 

sets.alignMethod = 1;
sets.edgeDetectionSettings.method = 'Otsu';

%% START CALCULATION
% load kymograph data
import Core.load_kymo_data;
[kymoStructs,barGen] = load_kymo_data(sets);

% save([sets.dirName, 'bars.mat'],'barGen','kymoStructs','sets');
save(['bars.mat'],'barGen','kymoStructs','sets');

% figure,tiledlayout(ceil(sqrt(length(kymoStructs))),ceil(length(kymoStructs)/sqrt(length(kymoStructs))),'TileSpacing','none','Padding','none')
% for i=1:length(kymoStructs)
%     nexttile;        imshowpair(imresize(kymoStructs{i}.alignedMask,[200 500]),imresize(kymoStructs{i}.alignedKymo,[200 500]), 'ColorChannels','red-cyan'  );    title(num2str(i));
% % imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  )
% end

% load(thryFile);
sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);

% extract from name
sets.theory.nmbp = sets.nmbp;

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct,sets );

% convert thry to correct nm/px

%
tic
sets.w = 300;
sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
globalov = 0;

import Core.extract_species_name; % find e-coli
[speciesLevel, idc] = extract_species_name(theoryStruct);
if globalov

% compare theory to experiment
import CBT.Hca.Core.Comparison.compare_distance;
[rezMax,bestBarStretch,bestLength] = compare_distance(barGen,theoryStruct, sets, [] );


%% Selected seq
% barIx = 1;
% toc
% allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMax);
% 
% [a,sortedid] = sort(allCCs,'desc','MissingPlacement','last');
% maxs = a(1);
% cdiff = 0.05;
% selectedRef = sortedid(find(a>a(1)-cdiff))

%% Local - how many significant matches



import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allNums,refNums,signMatch, fp,positives] = discrim_true_positives(rezMax, speciesLevel, idc);
theoryStruct([cell2mat(refNums(signMatch))]).name


theoryStruct([cell2mat(refNums(5))]).name

% also export info about disc species
import Core.export_coefs;
export_coefs(theoryStruct,rezMax,bestBarStretch,barGen,[sets.dirName, '_PCC_']);
end
% coefs to rezmax 

% import CBT.Hca.Export.export_cc_vals_table;
% [T] = export_cc_vals_table( theoryStruct, comparisonStructAll, barcodeGenC,matDirpath);


%%

% theoryStruct([refNumsMP{5}]).name;
% windowWidths = 400:100:600;
sets.comparisonMethod = 'mpnan';


import CBT.Hca.Core.Comparison.compare_distance;

for wIdx = 1:length(windowWidths)
    sets.w = windowWidths(wIdx);
    % only local lengths for which all length re-scaled versions passes the
    % threshold
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),barGen)*sets.theory.stretchFactors(end) >= sets.w);

    if sets.w == 0
        sets.comparisonMethod = 'mass_pcc';
    else
        sets.comparisonMethod = 'mpnan';
    end
    % assign standard scores
%     rezMaxMP = rezMax;
%     bestBarStretchMP = bestBarStretch;
%     bestLengthMP = bestLength;

    [rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen(passingThreshBars),theoryStruct, sets, [] );

    import Core.discrim_true_positives;
    [truePositivesMP{wIdx},discSpeciesMP{wIdx},discAllMP{wIdx},allNumsMP{wIdx},refNumsMP{wIdx},signMatchMP{wIdx}] =...
        discrim_true_positives(rezMaxMP,speciesLevel,idc);

%     sum(discSpecies(passingThreshBars)==0)
%     truePositivesMP{wIdx}

    %
    import Core.export_coefs;
    export_coefs(theoryStruct,rezMaxMP,bestBarStretchMP,barGen(passingThreshBars),[sets.dirName, '_MP_w=',num2str(sets.w),'_']);
    save([sets.dirName, num2str(sets.w),'_rez.mat'],'rezMaxMP','passingThreshBars','sets');

%     discSpecies(passingThreshBars)==0
% discSpeciesMP{wIdx}==0
% [truePositives,discSpecies,discAll,allNums,refNums,signMatch] = discrim_true_positives(rezMax,barGen,speciesLevel);
end


mpcalc = 0;
if mpcalc
% save output PCC
%% MPNAN
tic
sets.w = 500;
sets.comparisonMethod = 'mpnan';
notpassingThresh = find(cellfun(@(x) sum(x.rawBitmask),barGen)< sets.w);

import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen,theoryStruct, sets, [] );
%%%
toc
import Core.disc_locs;
[refNumsMP,allNumsMP] = disc_locs(rezMaxMP,barGen)
discAllMP = cellfun(@(x) ismember(x,allSpecies),refNumsMP,'UniformOutput',false)
discSpeciesMP = cellfun(@(x) sum(ismember(x,allSpecies)==0),refNumsMP,'UniformOutput',true)
sum(discSpeciesMP==0)

discSpeciesMP(notpassingThresh) = discSpecies(notpassingThresh)
sum(discSpeciesMP==0)



% refNumsMP(notpassingThresh)

theoryStruct([refNumsMP{5}]).name;

%
import Core.export_coefs;
export_coefs(theoryStruct,rezMaxMP,bestBarStretchMP,barGen,[sets.dirName, '_MP_w=',num2str(sets.w),'_']);

%% also MP to full for each

%% ASSEMBLY
import Core.run_local_assembly




%% Local - how many significant matches
refNumsMP = cell(1,length(barGen));
for barIx = 1:length(barGen);
%     toc
    allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMaxMP);

% figure,plot(allCCs)
% xlabel('Theory bar')
% ylabel('Score')

    [a,sortedid] = sort(allCCs,'desc','MissingPlacement','last');
    maxs = a(1);
    cdiff = 0.05;
    selectedRef = sortedid(find(a>a(1)-cdiff));
    refNumsMP{barIx} = selectedRef;
end


allNumsMP = cellfun(@(x) length(x),refNumsMP)
signMatchMP = find(allNumsMP ==1)
% refidsMP =  cell2mat(refNumsMP(signMatchMP))
% theoryStruct(refNumsMP{signMatchMP})
theoryStruct([refNumsMP{signMatchMP}]).name

[allNumsMP;allNums]

end

% quick_visual_plot(16,9242,barGen,rezMax,bestBarStretch,theoryStruct)

%  super_quick_plot(16,barGen,comparisonStruct,theoryStruct)
% sigmatches = find(allNums ==1)
% for i=1:length(sigmatches)
%     quick_visual_plot(sigmatches(i),9242,barGen,rezMax,bestBarStretch,theoryStruct)
% end

% cell2mat(refNums(sigmatches))
% refNums(signMatch)
% theoryStruct([cell2mat(refNums(signMatch))]).name;



% local_alignment_gui(sets)



% files = dir('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/Ecoli_DA65808/RawKymographs 1908XX - 110nm - 0.267nm_bp/*.tif');




end

