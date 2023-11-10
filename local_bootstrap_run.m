function [scores,pccScore] = local_bootstrap_run(barGenRun,rM,bnames,theoryStruct ,mpval,speciesLevel,idc)
% local_bootstrap_run
%   Args:
%       barGenRun, rM, bnames, theoryStruct , mpval, speciesLevel, idc   
%   Returns:
%       scores - mean( pccScore{m}) - 3*std(pccScore{m});
N = length(rM); % num different lengths

scores = nan(N,7);
pccScore = cell(1,N);

for m = 1:N % which length to analyse
    % Step 1 : find the best scoring theory for each barcode  
    import Core.disc_locs;
    [refNums, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{m});

%     % convert bestCoefs to p-values/zscores
    lengthsbar = sum(barGenRun{1}.rawBitmask);
%     
%     % pvalue params. should change based on psf
%     ccthresh = 0.3;
%     pval = nan(1,length(barcodeGen));
%     for i=1:length(barcodeGen)
% 
%         if comparisonStruct{i}.maxcoef(1) > ccthresh
%             pval(i) = 1-beta_ev_cdf(comparisonStruct{i}.maxcoef(1), nuF*lengthsbar(i), 1, 2*max(lengthsbar(i),lengthsThry(i)),1);
%         end
%     end



    % calcs ref and allnums again..
    import Core.discrim_true_positives;
        [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positives] = ...
        discrim_true_positives(rM{m}, speciesLevel, idc);

    % intermediate step, check if this barcode is within the output struc
    % (depeends on its length)
    [barsPassThresh,locb] = ismember(cellfun(@(x) matlab.lang.makeValidName(strrep(x.name,'.tif','')),barGenRun,'un',false),bnames{m});

    if barsPassThresh~=0 && sum(bestCoefs{locb})~=0

        % p-value for highest coef
        import Zeromodel.beta_ev_cdf;
        nuF = 0.05; % should be tuned based on synthetic data? Depends a bit on pdf
%         pvalAll = zeros(1,length(bestCoefs{locb}));
%         for i=1:length(bestCoefs{locb})
%             pvalAll(i) = 1-beta_ev_cdf(bestCoefs{locb}(i), nuF*rM{m}{refNums{locb}(i)}{locb}.lengthMatch, 1, 2*max(lengthsbar, theoryStruct(refNums{locb}(i)).length),1);
%         end
        pval1 = 1-beta_ev_cdf(bestCoefs{locb}(1), nuF*rM{m}{refNums{locb}(1)}{locb}.lengthMatch, 1, 2*max(lengthsbar, theoryStruct(refNums{locb}(1)).length),1);

        
        % Step 2 : recalculate (if needed) best alignment against theory. Useful to check
        % that the alignment was correct. Not needed if rezults is saved as mat
        % file
        theoryStructSel = theoryStruct(refNums{locb}(1)); % against the same theory, but we don't know which is correct! so we calcualte it for each length separately
        % barGenRun = barGen(1);
        % w = [200:50:sum(barGenRun{1}.rawBitmask)]; % in practice could run all
        [rezOutRecalc] = local_alignment_assembly(theoryStructSel, barGenRun,mpval(m));
    
        %


        % Step 3 : block bootstrapping
        [pccScore{m}, bar1, bar2] = block_bootstrapping(barGenRun, theoryStructSel,rezOutRecalc{1},1, 2,mpval(m));
%         figure,plot(zscore(bar1,1));
%         hold on
%         plot(zscore(bar2,1))
% %     

    
        
        pval = arrayfun(@(x) 1-beta_ev_cdf(x, nuF*rM{m}{refNums{locb}(1)}{locb}.lengthMatch, 1, 2*max(lengthsbar, theoryStruct(refNums{locb}(1)).length),0),pccScore{m});
        if min(pval) > 0.1 % consider only the case when pval > 0.01
            match_score = nan;
        else
            pval(pval>0.9999999999999) =0.9999999999999;
            pval(pval<10^(-20)) = 10^(-20);
            match_score = -norminv(pval);
        end

        scores(m,:) = [mean( pccScore{m}) 3*std(pccScore{m}) length(refNums{locb}) positives(locb) -norminv(pval1) mean(match_score) std(match_score)];
    else
    end
end

% figure,histogram(match_score)
% xlabel('match score (-norminv(pval,0,1)')
% ylabel('match score histogram')


end

