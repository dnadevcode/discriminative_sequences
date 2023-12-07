function [] = new_scores_pipeline(pathMain,dirName, depth, ix)
% pipeline for new scores (stouffer and bootstrapped)

%   Args:
%       pathMain - main path for reps
%       dirName - directory name with result files
%       depth - whether folders has subfolder structure, default 0


if isempty(pathMain) || nargin < 1
    pathMain = '/proj/snic2022-5-384/users/x_albdv/reps/';% reps.
end


if nargin < 2 || isempty(dirName)
    dirName = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/new_scores/';
end

if nargin < 3
    depth = 0;
end

addpath(genpath(dirName));
[~,~] =mkdir('output');


addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))


thryFiles = [dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat'),dir('/export/scratch/albertas/data_temp/bargrouping/New Ref-Theory files May 2022/*.mat')];


dirStruct = dir(dirName);
dirStruct(~[dirStruct.isdir]) = [];  %remove non-directories
dirStruct(ismember( {dirStruct.name}, {'.', '..'})) = [];  %remove . and ..




% for ix = 1;
%% - 
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

% this contains a single folder (containing kymographs). we can extract
% nm/bp of these
subDirKymos = dir(fullfile(sets.dirName ));
subDirKymos(ismember( {subDirKymos.name}, {'.', '..'})) = [];  %remove . and ..
subDirKymos(find(~cell2mat({subDirKymos.isdir}))) = [];


sets.dirNameKymos = fullfile(subDirKymos(1).folder,subDirKymos(1).name);

spltName = strsplit(sets.dirNameKymos ,'_');
spltName2 = strsplit(spltName{end},'nm');
spltName3 = strsplit(spltName{end-1},'nm');
sets.nmbp = str2double(spltName2{1});
nmpx = str2double(spltName3{1});


% load theory
thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,spltName{end-1}(1:5))),1:length(thryFiles)));
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);


% 
% curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/test/sample365/20220121_Sample365-st1_110nm_0.225nmPERbp/';

import Core.load_local_alignment_results_from_files;
[rM, bnames, mpval,thryNames,files] = load_local_alignment_results_from_files(sets.dirNameKymos ); 

%% if we do resampling, load theory
resampling = 1;

if resampling ==1
    import Core.load_theory_structure;
    [theoryStruct,sets2] = load_theory_structure(sets.nmbp,thryFileIdx);
end


[rM] = stouffer_zscore(rM,theoryStruct,mpval,0.07);

% load barcodes

fold = dir(fullfile(sets.dirName,'*bars.mat'));

barsfold = fullfile(fold(1).folder,fold(1).name);

bars=load(barsfold);


newFold = fullfile(strrep(dirName,'new_score','new_score_calc'),subDir(iy).name);
mkdir(newFold)
% export txts with stouffer data
import Core.export_coefs_stouffer;

for j=1:length(mpval)
    % slightly differed which barcodes to include
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen)*bars.sets.theory.stretchFactors(1) >= mpval(j));
    if length(passingThreshBars)~= length(rM{j}{1})
            passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen)*bars.sets.theory.stretchFactors(end) >= mpval(j));
    end
    if length(passingThreshBars)~= length(rM{j}{1})
        passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen) >=  mpval(j));
    end

    export_coefs_stouffer(theoryStruct,rM{j},bars.barGen(passingThreshBars),fullfile(newFold, ['stouffer_MP_w=',num2str(mpval(j)),'_']));
end

    




import Core.Default.read_default_sets;
hcaSets = read_default_sets('hcaalignmentsets.txt');
hcaSets.default.kymofolder{1} = 'kymo';
hcaSets.default.timestamp ='test';

import Core.identify_discriminative;
is_distinct = cell(1,length(rM));
numMatchingSpecies = cell(1,length(rM));
uniqueMatchSequences = cell(1,length(rM));
refNums = cell(1,length(rM));
refNumBad = cell(1,length(rM));
hcsets = hcaSets.default;

parfor m=1:length(rM)
    [is_distinct{m},numMatchingSpecies{m},uniqueMatchSequences{m},refNums{m},refNumBad{m}] =...
    identify_discriminative(hcsets, [], rM{m}, thryNames{m});
end



scores = cell(1,length(bars.barGen));

bg = bars.barGen;
for ii=1:length(bars.barGen)
    [scores{ii},~] = local_bootstrap_run( bg(ii),rM,bnames,theoryStruct ,is_distinct,numMatchingSpecies,uniqueMatchSequences,refNums, mpval);

end


[rM] = stouffer_zscore_bootstrap(rM,scores);



newFold = fullfile(strrep(dirName,'new_score','new_score_resamp'),subDir(iy).name);
mkdir(newFold)
% export txts with stouffer data
for j=1:length(mpval)
    % slightly differed which barcodes to include
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen)*bars.sets.theory.stretchFactors(1) >= mpval(j));
    if length(passingThreshBars)~= length(rM{j}{1})
            passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen)*bars.sets.theory.stretchFactors(end) >= mpval(j));
    end
    if length(passingThreshBars)~= length(rM{j}{1})
        passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),bars.barGen) >=  mpval(j));
    end

    import Core.export_coefs_stouffer_norm;                   
    export_coefs_stouffer_norm(theoryStruct,rM{j},bars.barGen(passingThreshBars),fullfile(newFold, ['stoufferresamp_MP_w=',num2str(mpval(j)),'_']));
end



% for m=1:length(rM)
%     for i=1:length(rM{m})
% 
%     end
% end


% baridxs = 1:length(rezMaxMP{i});
% baridxs = 1;
% for i =1:length(rezMaxMP)
%     i
%     for j=baridxs
%         % calc p-vals
%         if rezMaxMP{i}{j}.maxcoef(1) > minCC 
%             rezMaxMP{i}{j}.stoufferScoreNorm = (scores{j}(5) - rezMaxMP{i}{j}.stoufferScore)/scores{j}(6);
%         else
%             rezMaxMP{i}{j}.stoufferScoreNorm = nan;
%         end
% %          rezMaxMP{i}{j}.length = lengthMatch;
%     end
% end

%%

% import Core.Default.read_default_sets;
% hcaSets = read_default_sets('hcaalignmentsets.txt');
% hcaSets.default.kymofolder{1} = 'kymo';
% hcaSets.default.timestamp ='test';
% % 
% import Core.identify_discriminative;
% [is_distinct,numMatchingSpecies,uniqueMatchSequences,refNums,refNumBad] =...
%     identify_discriminative(hcaSets.default,barGen,rezMaxMP, thryNames{1});
%     


end

