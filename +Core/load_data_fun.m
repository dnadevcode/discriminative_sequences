function [kymoStructs,barGen,sets] = load_data_fun(dirName,ix,depth)


% EC17 = DA65780
% EC15 = DA65783
% EC6 = DA65808
% EC24 = DA65788
% EC25 = DA65793

% if nargin < 4
%     windowWidths = 400:100:600;
% end

if nargin < 3
    depth = 1;
end
% folds
% directory name
% thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat');

if nargin < 1
    dirName = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/';
end
addpath(genpath(dirName));
[~,~] =mkdir('output');
if nargin < 2
    ix = 1;
end

dirStruct = dir(dirName);
dirStruct(~[dirStruct.isdir]) = [];  %remove non-directories
dirStruct(ismember( {dirStruct.name}, {'.', '..'})) = [];  %remove . and ..


if depth == 1
    subDir = dir(fullfile(dirStruct(ix).folder,dirStruct(ix).name));
    subDir(ismember( {subDir.name}, {'.', '..'})) = [];  %remove . and ..
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
% thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,spltName{end-1}(1:5))),1:length(thryFiles)));
% sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);

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

end


