% resampling script

pathMain = 'C:\Users\Lenovo\git\dnadevcode\'; % path to reps.

foldToRun = '';

addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))

addpath('/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/')
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/'))


%dirName = 'C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment';
dirName = '/export/scratch/albertas/data_temp/bargrouping/local_results_from_tetralith/'; % work pc%
import Core.load_data_fun;
[kymoStructs,barGen,sets] = load_data_fun(dirName,5);

curDirKymos = sets.dirName; % current directory with kymographs

nmbp = sets.nmbp;

%% theory loading

import Core.load_theory_structure;
thryFileIdx = 1; % todo: pass directly the theory file here
[theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);


import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStruct);
%

%% Run alignment, if it has not been run yet

barGenRun = barGen(1);
% w = [];
% [rezMax, bestBarStretch, bestLength, rezOut] = local_alignment_assembly(theoryStruct, barGenRun,w);

%  
% idx = 1;
% import Core.discrim_true_positives;
% [truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
% discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);
% 
% {theoryStruct([refNums{1}]).name}'


%%

% foldCalc = 'C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment\EC6_Ecoli_DA65808';
% foldCalc = 'C:\Users\Lenovo\postdoc\DATA\bargrouping\local\local_alignment\local_alignment\EC15_Ecoli_DA65783';

% load alignment result to result struct if it was saves as text files
import Core.load_local_alignment_results_from_files;
[rM, bnames, mpval] = load_local_alignment_results_from_files(curDirKymos ); 


%% now run re-sampling procedure
% [scores,pccScore] = local_bootstrap_run( barGen(1),rM,bnames,theoryStruct ,mpval,speciesLevel,idc);
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

scores = cell(1,length(barGen));

for ii=1:length(barGen)
    [scores{ii},pccScore] = local_bootstrap_run( barGen(ii),rM,bnames,theoryStruct ,mpval,speciesLevel,idc);
    scores{ii}
    
    % output:
    import Core.export_coefs_resampling;
    T = export_coefs_resampling(scores{ii}, barGen(ii), mpval, [curDirKymos, '_resampling_table'],timestamp);
end

%% same individually

m = 10; % which length to analyse
% Step 1 : find the best scoring theory for each barcode  
import Core.disc_locs;
[refNums, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{m});

% intermediate step, check if this barcode is within the output struc
% (depeends on its length)
[barsPassThresh,locb] = ismember(cellfun(@(x) matlab.lang.makeValidName(strrep(x.name,'.tif','')),barGenRun,'un',false),bnames{m});

% Step 2 : recalculate (if needed) best alignment against theory. Useful to check
% that the alignment was correct. Not needed if rezults is saved as mat
% file

theoryStructSel = theoryStruct(refNums{locb}(1)); % against the same theory, but we don't know which is correct! so we calcualte it for each length separately
% barGenRun = barGen(1);
% w = [200:50:sum(barGenRun{1}.rawBitmask)]; % in practice could run all
[rezOutRecalc] = local_alignment_assembly(theoryStructSel, barGenRun,mpval(m));

%     super_quick_plot(1,barGenRun,rezOutRecalc{1},theoryStructSel)
%% super_quick_plot_mp 
% import Core.plot_match_simple;
% % [f] = plot_match_simple(barStruct, oS,curSink,curSource);
% [f] = plot_match_simple([barGenRun(selRef) barcodeGenT],rezOut{2}, 1, refNums{selRef}(idx)+1);
% 

% Step 3 : block bootstrapping
pccScore = block_bootstrapping(barGenRun, theoryStructSel,rezOutRecalc{1},1, 2,mpval(m));

stdVals(k-1) = mean( pccScore) - 3*std(pccScore);

%%



figure,plot(mpval,scores(:,1)-3*(scores(:,2)),'redx')
xlabel('Overlap window width');

ylabel('$CC_{max}$ -3$\sigma_{bootstrapped}$','Interpreter','latex')

%% 

% local_bootstrap_run(barGenRun,rM,bnames,theoryStruct,mpval )
% 
% stdVals = zeros(1,length(rezOut)-1);
% for k=2:length(rezOut)
%  
%     barcodeGenTStruct = struct('rawBarcode',theoryStructNew(1).rawBarcode,'rawBitmask',true(1,length(theoryStructNew(1).rawBarcode)));
% 
%     pccScore = block_bootstrapping([barGenRunStruct(selRef) barcodeGenTStruct(1)],rezOut1{k},1, 2);
%     stdVals(k-1) = mean( pccScore) - 3*std(pccScore);
% end




import Core.discrim_true_positives;
pos = zeros(1,length(files));
for i=1:length(files)
    [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positivesMP] = ...
        discrim_true_positives(rM{i}, speciesLevel, idc);
    pos(i) = length(cell2mat(refNumsMP));
end


%

theoryStructSel = theoryStruct(refNums{1}(1));
barGenRun = bgAll(1);
w = [200:50:sum(barGenRun{1}.rawBitmask)];
[rezMax, bestBarStretch, bestLength, rezOut1] = local_alignment_assembly(theoryStructSel, barGenRun,w);

theoryStructSel2 = theoryStruct(refNums{1}(2:end));
[rezMax2, bestBarStretch2, bestLength2, rezOut2] = local_alignment_assembly(theoryStructSel2, barGenRun,w);

%
 barGenRunStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barGenRun,'un',false);...
        cellfun(@(x) x.rawBitmask,barGenRun,'un',false)]',{'rawBarcode','rawBitmask'},2);
% import Core.plot_match_pcc;
% [sortedValsIs, sortedIdsIs,pscoresIs] = sorted_scores(oSIslands);




%% Other things




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



