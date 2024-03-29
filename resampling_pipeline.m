% resampling pipeline

% We're testing all the different possible methods to get a standard error
% on the local alignment scores

% Should run this for global and local alignment

addpath(genpath('/home/avesta/albertas/reps/lldev'))
addpath(genpath('/home/avesta/albertas/reps/hca'))
addpath(genpath('/home/avesta/albertas/reps/bargroupingprototype'))
addpath(genpath('/home/avesta/albertas/reps/discriminative_sequences'))

addpath('/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/')
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/'))




% data loading
import Core.load_chrom_data;
[bgAll, bG, kymoStructs] = load_chrom_data('/export/scratch/albertas/data_temp/bargrouping_selected/ecoli_1/sample_2/');
fastaFiles(1:length(bG)) = 1;

nmbp = 0.25;

%% theory loading

import Core.load_theory_structure;
thryFileIdx = 1; % todo: pass directly the theory file here
[theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);

%% optional: keep only e-coli

%%
%% todo: also add current theoretical barcode to theoryStruct
import Thry.gen_theoretical;

nmPerPx = 110;
fastas = {'018_final_polish.fasta','DA32087.fasta'};

fastaFiles = 1:2;

% if for all nm/bp values
nmbpvals = 0.25; % should use mpMax to get this "fully" correct
[theoryStructNew,~,barcodeGenT] = gen_theoretical(fastas(fastaFiles),nmbpvals,0,nmPerPx); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly

%% add new theory to the end
theoryStruct = [theoryStruct; theoryStructNew];

% theoryStruct = theoryStructNew
%% compare and get discriminative stuff

barGenRun = bgAll(1:10);%(1);
w = [];
[rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStruct, barGenRun,w);


import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStruct);
% 
idx = 1;
import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);

{theoryStruct([refNums{4}]).name}'

%% visualize the best result
selRef = 5;
idx = 1;
idx1 = selRef;
quick_visual_plot(1,refNums{selRef}(idx),barGenRun,rezMax,bestBarStretch,theoryStruct)
super_quick_plot(5,barGenRun,rezOut{1},theoryStruct)
%% super_quick_plot_mp
import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple([barGenRun(selRef) barcodeGenT],rezOut{2}, 1, refNums{selRef}(idx)+1);

%% first
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

stdVals = zeros(1,length(rezOut)-1);
for k=2:length(rezOut)
 
    barcodeGenTStruct = struct('rawBarcode',theoryStructNew(1).rawBarcode,'rawBitmask',true(1,length(theoryStructNew(1).rawBarcode)));

    pccScore = block_bootstrapping([barGenRunStruct(selRef) barcodeGenTStruct(1)],rezOut1{k},1, 2);
    stdVals(k-1) = mean( pccScore) - 3*std(pccScore);
end

% figure,plot(w,stdVals)
% xlabel('window width (w)')
% ylabel('mu-3sigma best barcode bootstrapping')


%second best
stdVals2 = zeros(1,length(rezOut)-1);
for k=2:length(rezOut)
    [~,~,~,~,refNums2,~] =...
    discrim_true_positives(rezOut2{k}.rezMax,speciesLevel,idc);
    barcodeGenTStruct = struct('rawBarcode',theoryStructSel2(refNums2{1}(1)).rawBarcode,'rawBitmask',true(1,length(theoryStructSel2(refNums2{1}(1)).rawBarcode)));

    barGenRunStruct(refNums2{1}(1)+1) = barcodeGenTStruct;
%     temps = [theoryStructSel2.rawBarcode,theoryStructSel2.rawBitmask];
%     barcodeGenTStruct = struct('rawBarcode',theoryStructSel2.rawBarcode);

    pccScore2 = block_bootstrapping(barGenRunStruct,rezOut2{k},1, refNums2{1}(1)+1);
    stdVals2(k-1) = mean( pccScore2) + 3*std(pccScore2);
end

% figure,plot(w,stdVals2)
% xlabel('window width (w)')
% ylabel('mu+3sigma 2nd best bootstrapping')



%second best
numElts = zeros(1,length(rezOut)-1);
for k=2:length(rezOut)
    [~,~,~,~,refNumsAll,~] =...
    discrim_true_positives([rezOut1{k}.rezMax rezOut2{k}.rezMax],speciesLevel,idc);
    numElts(k-1) = length(refNumsAll{1});
%     pccScore2 = block_bootstrapping([barGenRun(1) barcodeGenTStruct],rezOut{k},1, refNums2{1}(1)+1);
%     stdVals2(k-1) = mean( pccScore2) + 3*std(pccScore2);
end


score = stdVals - stdVals2;
figure,tiledlayout(2,1);
nexttile
plot(w,score)
hold on
xlabel('window width (w)')
ylabel('score')
title('$score = (\mu-3\sigma)_{best} -(\mu+3\sigma)_{2ndbest}$','Interpreter','latex')
nexttile
plot(w,numElts)
xlabel('window width (w)')
title('Number of matching theories within 0.05')
%%
tic
% barGenRun = bgAll;
% sets.theory.stretchFactors = 0.9:0.025:1.1; %as per 
% sets.w = 500;
% sets.comparisonMethod = 'mass_pcc'; % mass_pcc
% sets.genConsensus = 0;
% sets.filterSettings.filter = 0;
% sets.dirName
% % barGen = bgAll(1:200);
% % compare theory to experiment
% import CBT.Hca.Core.Comparison.compare_distance;
% [rezMaxNew,bestBarStretchNew,bestLengthNew] = compare_distance(barGenRun,theoryStructNew, sets, [] );

maxCC = [cellfun(@(x) x.maxcoef(1),rezMax{2})];
% [mean([cellfun(@(x) x.maxcoef(1),rezMax{1})]) std([cellfun(@(x) x.maxcoef(1),rezMax{1})])]
maxCC = arrayfun(@(y) rezMax{end-1}{y}.maxcoef(1),1:length(refNums),'UniformOutput',true)


theoryStruct(end+1) = theoryStructNew;
%%

% for best
maxCoefsBest = arrayfun(@(y) max(arrayfun(@(x) max(rezMax{x}{y}.maxcoef(1)),refNums{y}(:))),1:length(refNums),'UniformOutput',true)

% difference between correct and second
maxCC(1:length(maxCoefsBest)) - maxCoefsBest

%%

rezMax = rezOut{1}.rezMax;
bestBarStretch = rezOut{1}.bestBarStretch;
% rezMax = rezOut{1}.rezMax;

selRef=1;

%     maxCoefs = arrayfun(@(x) rezMax{x}{selRef}.maxcoef(1),refNums{selRef}(:));

maxCoefs = arrayfun(@(x) rezMax{x}{selRef}.maxcoef(1),1:length(rezMax));

[sortMax,sortMaxId] = sort(maxCoefs,'desc');

    idx = 1;
    idx1 = selRef;
    quick_visual_plot(1,refNums{selRef}(idx),barGenRun,rezMax,bestBarStretch,theoryStruct)
    super_quick_plot(1,barGenRun,rezOut{1},theoryStruct)

    % time-frame bootstrapping
    timeframe_bootstrapping

    %%
thrIdx = refNums{selRef}(idx);%length(theoryStruct)-1; %refNums{1}(idx);
    
    lenBarTested = length(barGenRun{idx1}.rawBarcode);
    bar1 = interp1(barGenRun{idx1}.rawBarcode, linspace(1,lenBarTested,lenBarTested*bestBarStretch{thrIdx}(idx1)));
    barB = barGenRun{idx1}.rawBitmask(round(linspace(1,lenBarTested,lenBarTested*bestBarStretch{thrIdx}(idx1))));
    
    
    % bar1 = imresize(barGen{idx1}.rawBarcode(barGen{idx1}.rawBitmask),'Scale' ,[1 bestBarStretch{thrIdx}(idx1)]) ;
    if rezMax{thrIdx}{idx1}.or(1) == 2
        bar1 = fliplr(bar1);
    end
    bar1 = bar1(barB);
    
    {theoryStruct([refNums{selRef}]).name}'
    
    thr = theoryStruct(thrIdx).theoryBarcode;
    
    pos = find(barGenRun{idx1}.rawBitmask==1,1,'first');
    bar2 = thr(rezMax{thrIdx}{idx1}.pos+pos-1:rezMax{thrIdx}{idx1}.pos+length(bar1)-1+pos-1);
    % 
    zscore(bar1(:)',1)*zscore(bar2(:),1)/length(bar1(:))
    
    figure,plot(zscore(bar1,1));
    hold on
    plot(zscore(bar2))
    
    % now split into ~50px windows
    wWid = 20;
    numWindows = floor(length(bar1)/wWid);
    
    R = reshape(1:floor(length(bar1)/numWindows)*numWindows,floor(length(bar1)/numWindows),[]);
%     s = RandStream('mlfg6331_64'); 
    
    NN = 1000;
    pccScore = zeros(1,NN);
    for i=1:NN
        y = datasample(s,1:numWindows,numWindows,'Replace',true);
        b1 = bar1(R(:,y));
        b2 = bar2(R(:,y));
        pccScore(i) = zscore(b1(:)')*zscore(b2(:))/length(b1(:));
    
    end
    % 
    % figure,plot(zscore(b1(:)));
    % hold on
    % plot(zscore(b2(:)))
    
    mean(pccScore)
    [3*std(pccScore) 3*std(sortMax(2:end))]
    
    figure,histogram(pccScore)
    xlabel('CC')


    figure,plot(zscore(bar1,1));
    hold on
    plot(zscore(bar2))
    
%%
    kymo =barGenRun{5}.alignedKymo;
    % figure,imagesc(kymo)
    meankymo = nanmean(kymo);
    stdkymo = nanstd(kymo);

    figure,
    tiledlayout(1,2)
    nexttile([1 2])
%     errorbar(    mean(kymo),stdkymo)

% 
stdkymo = stdkymo(~isnan(meankymo));
meankymo = meankymo(~isnan(meankymo));
x = 1:length(meankymo);

xconf = [1:length(meankymo) length(meankymo):-1:1] ;  
    y =  meankymo;
yconf = [y+3*stdkymo y(end:-1:1)-3*stdkymo(end:-1:1)];
% 
% figure
p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';           
% 
hold on
plot(x,y,'r-')
hold off



    %% Same but local alignment
    w = 300;


    %%

import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStructSel);
% 
idx = 10;
import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allSpecies,refNumsSel,signMatch] =...
discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);

{theoryStructSel([refNumsSel{1}]).name}'