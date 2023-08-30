
% Should run this for global and local alignment

addpath(genpath('/home/etlar/albertas/reps/lldev'))
addpath(genpath('/home/etlar/albertas/reps/hca'))
addpath(genpath('/home/etlar/albertas/reps/bargroupingprototype'))
addpath(genpath('/home/etlar/albertas/reps/discriminative_sequences'))

addpath('/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/')
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/'))



% data loading
import Core.load_chrom_data;
[bgAll, bG, kymoStructs] = load_chrom_data('/export/scratch/albertas/data_temp/hca_data/');
% [bgAll, bG, kymoStructs] = load_chrom_data('/export/scratch/albertas/data_temp/tests/');

fastaFiles(1:length(bG)) = 1;

% nmbp = 0.225;

lens = cellfun(@(x) length(x.rawBarcode),bgAll)
nmbp = 0.22;
bgsel = bgAll(find(lens>500));

%% theory loading

import Core.load_theory_structure;
thryFileIdx = 1; % todo: pass directly the theory file here
[theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);

%%

barGenRun = bgsel(3);

w = [];
[rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStruct, barGenRun,w);


w = [200:50:sum(barGenRun{1}.rawBitmask)];
[rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStruct, barGenRun,w);




import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStruct,{'Streptococcus'});
% 
idx = 1;
import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);

{theoryStruct([refNums{1}]).name}'

%%
barGenRun = bgAll(1);
w = [];
[rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStruct, barGenRun,w);


selRef = 1;
idx = 1;
idx1 = selRef;
quick_visual_plot(1,refNums{selRef}(idx),barGenRun,rezMax,bestBarStretch,theoryStruct)
super_quick_plot(1,barGenRun,rezOut{1},theoryStruct)
%% super_quick_plot_mp
import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple([barGenRun(selRef) barcodeGenT],rezOut{2}, 1, refNums{selRef}(idx)+1);




%% only single species
theoryStructSel = theoryStruct(find(speciesLevel));
barGenRun = bgAll(1);
w = [];
[rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStructSel, barGenRun,w);

figure,plot(cellfun(@(x) x{1}.maxcoef(1),rezMax))
w = [200:50:sum(barGenRun{1}.rawBitmask)];
[rezMax, bestBarStretch, bestLength, rezOut1] = local_alignment_assembly(theoryStructSel, barGenRun,w);

% check nm/bp ratio



%% nm bp
% Evaluation figure - experiments from individual days
thryFiles = dir('/export/scratch/albertas/data_temp/bargrouping/New Ref-Theory files May 2022/*.mat');

thryFileIdx = 1;
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);


sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStructA = load_theory(sets);


%% from here independent code
theoryStructInd =theoryStructA((find(speciesLevel)));


sets.comparisonMethod = 'mass_pcc';
sF = 0.9:0.01:1.1;
minLen = 300;


%% from here independent code
% theoryStructInd =theoryStruct([refNums{1}]);
import CBT.Hca.Core.Analysis.convert_nm_ratio;
barGenRun = bgAll(1);

% psffac = 1; % scaling factor (in case PSF needs to be something else than 300nm)
numWorkers = 30; % num workers, in one node there are 30
% minLen = [150:50:3000]; % min to max % could take more points for more accurate..
sets.comparisonMethod = 'mass_pcc';
sF = 0.9:0.01:1.1;
%
% %     % check for different nmbp (in case estimated incorrectly
nmbpvals = 0.15:0.01:0.3; % should use mpMax to get this "fully" correct
import Thry.gen_theoretical;
scoef = zeros(1,length(nmbpvals));
for nmIdx =1:length(nmbpvals)
    tic
        theoryStructAll{nmIdx} = convert_nm_ratio(nmbpvals(nmIdx),theoryStructInd ,sets );

%                 [theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbpvals(nmIdx),0,nmPerPx); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
    [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(barGenRun,  theoryStructAll{nmIdx} ,sF,sets);
    
    allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
    allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
    
%         idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
%         coefDiffs = allCoefs-mpMaxLenBasedAll{nmIdx}(idxThresh);
    
    scoef(nmIdx) = allCoefs;
%         succesNmnp(nmIdx) = length(barsPassThresh)/length(coefDiffs); % this could still have false positives.. ?
%         succesNmnp(nmIdx) = length(find(cellfun(@(x) x.maxcoef(1),comparisonStruct)-mpMaxLenBasedAll{nmIdx}(idxThresh) > 0.05))/length(bG{idxRun}); % this could still have false positives.. ?
    toc
end
%     sucRateStruct{idxRun}.succesNmnp = succesNmnp;
% 
    f = figure
    plot(nmbpvals,scoef);


    %%
    figure,histogram(cellfun(@(x) std(x.rawBarcode),bgAll)./cellfun(@(x) nanmean(std(x.alignedKymo)),bgAll))
