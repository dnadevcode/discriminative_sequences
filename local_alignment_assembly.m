function [rezMax,bestBarStretch,bestLength,discSpecies,refNums] = local_alignment_assembly(theoryStruct, barGen,w)
    % Tailored to work with kymographs from barcode assembly problem to
    % find which are discriminative, so we skip to bargen

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sets.dirName = 'output';

% 
% 
if nargin < 1
    % create sets
end

% import Core.load_theory_structure;
% thryFileIdx = 1; % todo: pass directly the theory file here
% [theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);


%  following "Strain-level bacterial typing directly from patient
% samples using optical DNA mapping"
sets.timeFramesNr = 20;
sets.theory.stretchFactors = 0.9:0.025:1.1; %as per 

sets.alignMethod = 1;
sets.edgeDetectionSettings.method = 'Otsu';


%
tic
sets.w = w;
sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
% barGen = bgAll(1:200);
% compare theory to experiment
import CBT.Hca.Core.Comparison.compare_distance;
[rezMax,bestBarStretch,bestLength] = compare_distance(barGen,theoryStruct, sets, [] );

% toc
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
[refNums,allNums] = disc_locs(rezMax)

signMatch = find(allNums ==1);
% refNums(signMatch)
% theoryStruct([cell2mat(refNums(signMatch))]).name;
% theoryStruct([cell2mat(refNums(:))]).name;

import Core.extract_species_name;
[speciesLevel] = extract_species_name(theoryStruct);

allSpecies = find(speciesLevel);

discAll = cellfun(@(x) ismember(x,allSpecies),refNums,'UniformOutput',false)

discSpecies = cellfun(@(x) sum(ismember(x,allSpecies)==0),refNums,'UniformOutput',true)
sum(discSpecies==0)
% 
% % also export info about disc species
% import Core.export_coefs;
% export_coefs(theoryStruct,rezMax,bestBarStretch,barGen,[sets.dirName, '_PCC_']);
% % import CBT.Hca.Export.export_cc_vals_table;
% % [T] = export_cc_vals_table( theoryStruct, comparisonStructAll, barcodeGenC,matDirpath);
% 


% theoryStruct([refNumsMP{5}]).name;
windowWidths = 400;%400:100:600;
sets.comparisonMethod = 'mpnan';


import CBT.Hca.Core.Comparison.compare_distance;

for wIdx = 1:length(windowWidths)
    sets.w = windowWidths(wIdx);
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),barGen) >= sets.w);

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

%     %
%     import Core.export_coefs;
%     export_coefs(theoryStruct,rezMaxMP,bestBarStretchMP,barGen(passingThreshBars),[sets.dirName, '_MP_w=',num2str(sets.w),'_']);

%     discSpecies(passingThreshBars)==0
% discSpeciesMP{wIdx}==0
% [truePositives,discSpecies,discAll,allNums,refNums,signMatch] = discrim_true_positives(rezMax,barGen,speciesLevel);
end


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

% cell2mat(refNums(sigmatches))
% refNums(signMatch)
% theoryStruct([cell2mat(refNums(signMatch))]).name;


end