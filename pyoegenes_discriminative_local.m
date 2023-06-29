% For pyogenes, run local comparison (from e-coli paper I)

% Evaluation figure - experiments from individual days
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/*.mat');

thryFileIdx = 1;
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);


sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);


%% from here independent code
theoryStructInd =theoryStruct(1129);
theoryStructT = convert_nm_ratio(sets.theory.nmbp,theoryStructInd ,sets );


import CBT.Hca.Core.Analysis.convert_nm_ratio;

% psffac = 1; % scaling factor (in case PSF needs to be something else than 300nm)
numWorkers = 30; % num workers, in one node there are 30
minLen = [150:50:3000]; % min to max % could take more points for more accurate..
sets.comparisonMethod = 'mass_pcc';
sF = 0.9:0.01:1.1;
nmPerPx = 110;

% fastas = {'018_final_polish.fasta','DA32087.fasta'};

% if for all nm/bp values
nmbpvals = 0.15:0.01:0.3; % should use mpMax to get this "fully" correct
import Thry.gen_theoretical;
theoryStructAll = cell(1,length(nmbpvals));
mpMaxLenBasedAll =  cell(1,length(nmbpvals));
for nmIdx =1:length(nmbpvals)
    nmIdx
    theoryStructAll{nmIdx} = convert_nm_ratio(sets.theory.nmbp,theoryStructInd ,sets );
    lenBp = theoryStructAll{nmIdx}.length* theoryStructAll{nmIdx}.pixelWidth_nm/theoryStructAll{nmIdx}.meanBpExt_nm;

%     [theoryStructAll{nmIdx},~,barcodeGenT] = gen_theoretical(fastas(fastaFiles(1)),nmbpvals(nmIdx),0,kymoStructs{1}{1}.nmpxidFold); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
    [~,mpMaxLenBasedAll{nmIdx},~,~,~] = bargrouping_minimum_length([],nmPerPx,nmbpvals(nmIdx),1,numWorkers,minLen, lenBp);
end

% f = figure;
% heatmap(minLen,nmbpvals,cell2mat(mpMaxLenBasedAll'))
% title('Heatmap for the local similarities threshold')
% print('FIGS/FigS1.eps','-depsc','-r300');

% xlabel('Overlap length','Interpreter','latex');
% ylabel('nm/bp extension factor','Interpreter','latex')

bG{1} = barGen;
nmbp = 0.169;

successRate = zeros(1,length(bG));
sucRateStruct = cell(1,length(bG));
comparisonStructC= cell(1,length(bG));
passthreshC= cell(1,length(bG));
mpMaxLenBasedC =  cell(1,length(bG));
for idxRun = 1:length(bG)
%     nmPerPx = kymoStructs{idxRun}{1}.nmpxidFold;
%     nmbp = kymoStructs{idxRun}{1}.nmBpidFold;
%     fastaFile = fastas(fastaFiles(idxRun));
    % can be calculated before for all nmbp
    %
    [MP,mpMaxLenBasedC{idxRun},theoryStructRev,MPI,~] = bargrouping_minimum_length([],nmPerPx,nmbp,1,numWorkers,minLen, theoryStructT.length*theoryStructT.pixelWidth_nm/theoryStructT.meanBpExt_nm);% no lenseq if based on fastafile

%     import Thry.gen_theoretical;
%     [theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);
tic
    import CBT.Hca.Core.Comparison.compare_distance;
    [rezMax2,bestBarStretch,bestLength] = compare_distance(bG{1},theoryStructT, sets, [] );
toc

% tic
    [comparisonStructC{idxRun},rezMax,bestBarStretch] = compare_to_t(bG{idxRun},{theoryStructT},sF,sets);
% toc
    allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStructC{idxRun});
    allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStructC{idxRun});
    
    idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
    coefDiffs = allCoefs-mpMaxLenBasedC{idxRun}(idxThresh);
    
    passthreshC{idxRun} = find(coefDiffs > 0.05);
    successRate(idxRun) = length(passthreshC{idxRun} )/length(coefDiffs); % this could still have false positives.. ?
    
% %     % check for different nmbp (in case estimated incorrectly
    nmbpvals = 0.15:0.01:0.3; % should use mpMax to get this "fully" correct
    import Thry.gen_theoretical;
    succesNmnp = zeros(1,length(nmbpvals));
    for nmIdx =1:length(nmbpvals)
        %         [theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbpvals(nmIdx),0,nmPerPx); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
        [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(bG{idxRun},{theoryStructAll{nmIdx}},sF,sets);
        
        allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
        allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
        
        idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
        coefDiffs = allCoefs-mpMaxLenBasedAll{nmIdx}(idxThresh);
        
        barsPassThresh = find(coefDiffs > 0.05);
        succesNmnp(nmIdx) = length(barsPassThresh)/length(coefDiffs); % this could still have false positives.. ?
%         succesNmnp(nmIdx) = length(find(cellfun(@(x) x.maxcoef(1),comparisonStruct)-mpMaxLenBasedAll{nmIdx}(idxThresh) > 0.05))/length(bG{idxRun}); % this could still have false positives.. ?
    end
    sucRateStruct{idxRun}.succesNmnp = succesNmnp;
% 
    f = figure
    plot(nmbpvals,succesNmnp);
%     xlabel('nm/bp','Interpreter','latex')
%     ylabel('succes rate','Interpreter','latex');
%     hold on
%     plot(nmbp,0,'redx')
%     lgd = legend({'Success rate','Estimated nm/bp'})
%     lgd.Location ='southoutside';
%     title('Success rate for different nm/bp extension factors','Interpreter','latex')
%     print('FIGS/FigS2.eps','-depsc','-r300');
% 
end
%%
idxRun = 9;
     f = figure
    plot(nmbpvals,   sucRateStruct{idxRun}.succesNmnp);
    xlabel('nm/bp','Interpreter','latex')
    ylabel('succes rate','Interpreter','latex');
    hold on
    plot(kymoStructs{idxRun}{1}.nmBpidFold,0,'redx')
    lgd = legend({'Success rate','Estimated nm/bp'})
    lgd.Location ='southoutside';
    title('Success rate for different nm/bp extension factors','Interpreter','latex')
%
%%

goodbars = bG{idxRun};
w = 300;
[rezMax,bestBarStretch,bestLength,discSpecies] = local_alignment_assembly(goodbars, nmbp,w);

% %%
% bpPx = 500;
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], comparisonStructC, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6
% % 
% bg = bgAll(barsPassThresh);
% [outConsensus, coverage, pval] = gen_reference_based_assembly(bg,comparisonStruct(barsPassThresh),theoryStruct,'test11',inf);
% 
% 
% % quick_visual_plot(46,1,bgAll,rezMax,bestBarStretch,theoryStruct)
% 
% [compStr,~,calcLengths] = compare_to_t_mp(bgAll,theoryStruct,sF,300); %todo: extend local to full
% [comparisonStruct{4}.pos(1) compStr{2}.pos(1)]
% 
% passthreshLen = find(ismember(calcLengths,barsPassThresh));
% % passthreshLen = barsPassThresh(calcLengths);
% 
% bpPx = 500;
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], compStr(passthreshLen), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6
% 
% % compStr{passthreshLen(1)}.pos = 4169;
% % super_quick_plot(1,bgAll(4),compStr(2),theoryStruct)
% 
% % allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% % allLengths =300;
% % 
% % idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
% % coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% % 
% % barsPassThresh = find(coefDiffs > 0);
% % 
% % bg = bgAll(barsPassThresh); % make gen_reference_based_assembly it work for MP
% [outConsensus, coverage, pval] = gen_reference_based_assembly(bgAll(passthreshLen),compStr(passthreshLen),theoryStruct,'test11',inf);
% 
