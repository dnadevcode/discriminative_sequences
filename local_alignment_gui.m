function [] = local_alignment_gui(sets)
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% 
% 
if nargin < 1
    % create sets
end

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
% load kymograph data
import Core.load_kymo_data;
[kymoStructs,barGen] = load_kymo_data(sets);


figure,tiledlayout(ceil(sqrt(length(kymoStructs))),ceil(length(kymoStructs)/sqrt(length(kymoStructs))),'TileSpacing','none','Padding','none')
for i=1:length(kymoStructs)
    nexttile;        imshowpair(imresize(kymoStructs{i}.alignedMask,[200 500]),imresize(kymoStructs{i}.alignedKymo,[200 500]), 'ColorChannels','red-cyan'  );    title(num2str(i));
% imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  )
end

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
import Core.disc_locs;
[refNums,allNums] = disc_locs(rezMax,barGen)

signMatch = find(allNums ==1)
% refNums(signMatch)
theoryStruct([cell2mat(refNums(signMatch))]).name;

import Core.extract_species_name;
[speciesLevel] = extract_species_name(theoryStruct);

allSpecies = find(speciesLevel);

discAll = cellfun(@(x) ismember(x,allSpecies),refNums,'UniformOutput',false)

discSpecies = cellfun(@(x) sum(ismember(x,allSpecies)==0),refNums,'UniformOutput',true)
sum(discSpecies==0)

% also export info about disc species
import Core.export_coefs;
export_coefs(theoryStruct,rezMax,bestBarStretch,barGen,[sets.dirName, '_PCC_']);
% import CBT.Hca.Export.export_cc_vals_table;
% [T] = export_cc_vals_table( theoryStruct, comparisonStructAll, barcodeGenC,matDirpath);
mpcalc = 0;
if mpcalc
% save output PCC
%% MPNAN
tic
sets.w = 200;
sets.comparisonMethod = 'mpnan';

import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen,theoryStruct, sets, [] );
%%%
toc
import Core.disc_locs;
[refNumsMP,allNumsMP] = disc_locs(rezMaxMP,barGen)
discAllMP = cellfun(@(x) ismember(x,allSpecies),refNumsMP,'UniformOutput',false)
discSpeciesMP = cellfun(@(x) sum(ismember(x,allSpecies)==0),refNumsMP,'UniformOutput',true)
sum(discSpeciesMP==0)

theoryStruct([refNumsMP{5}]).name

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

cell2mat(refNums(sigmatches))
% refNums(signMatch)
% theoryStruct([cell2mat(refNums(signMatch))]).name;


end