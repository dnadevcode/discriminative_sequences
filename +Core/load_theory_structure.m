function [theoryStruct,sets] = load_theory_structure(nmbp, thryFileIdx,selected)


sets.dirName = 'output';
sets.nmbp = nmbp;
thryFiles = [dir('/export/scratch/albertas/data_temp/bargrouping/New Ref-Theory files May 2022/*.mat'),dir('C:\Users\Lenovo\postdoc\DATA\bargrouping\*.mat')];
% thryFileIdx = 1;
% thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,spltName{end-1})),1:length(thryFiles)));
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);

% files = dir(fullfile(sets.dirName,'*.tif'));

% sets.kymosets.filenames = arrayfun(@(x) files(x).name,1:length(files),'un',false);
% sets.kymosets.kymofilefold = arrayfun(@(x) files(x).folder,1:length(files),'un',false);

% simdata = 0;
sets.output.matDirpath = 'output';
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 6;
sets.minLen = 250;
sets.genConsensus  = 0;


% load kymograph data
% import Core.load_kymo_data;
% [kymoStructs,barGen] = load_kymo_data(sets);

% 
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

if nargin >=3
    theoryStruct = theoryStruct(selected);
end

% extract from name
sets.theory.nmbp = sets.nmbp;

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct,sets );

% convert thry to correct nm/px


end

