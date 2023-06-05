
% EC17 = DA65780
% EC15 = DA65783
% EC6 = DA65808
% EC24 = DA65788
% EC25 = DA65793

% directory name
dirName = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/';
dirStruct = dir(dirName);
dirStruct(~[dirStruct.isdir]) = [];  %remove non-directories
dirStruct(ismember( {dirStruct.name}, {'.', '..'})) = [];  %remove . and ..

ix = 1;


subDir = dir(fullfile(dirStruct(ix).folder,dirStruct(ix).name));
subDir(ismember( {subDir.name}, {'.', '..'})) = [];  %remove . and ..

%
iy = 1; % most likely single run
sets.dirName = fullfile(subDir(iy).folder,subDir(iy).name);

spltName = strsplit(sets.dirName ,'_');
spltName2 = strsplit(spltName{end},'nm');
spltName3 = strsplit(spltName{end-1},'nm');
sets.nmbp = str2double(spltName2{1});
nmpx = str2double(spltName3{1});


% load theory
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat');
thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,spltName{end-1})),1:length(thryFiles)));
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);
local_alignment_gui(sets)



% files = dir('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/Ecoli_DA65808/RawKymographs 1908XX - 110nm - 0.267nm_bp/*.tif');
