function [scores,pccScore] = local_bootstrap_run(barGenRun,rM,bnames,theoryStruct ,is_distinct,numMatchingSpecies,uniqueMatchSequences,refNums, mpval)
% local_bootstrap_run
%   Args:
%       barGenRun, rM, bnames, theoryStruct , mpval, speciesLevel, idc   
%   Returns:
%       scores - mean( pccScore{m}) - 3*std(pccScore{m});
N = length(rM); % num different lengths

scores = nan(N,7);
pccScore = cell(1,N);

for m = 1:N % which length to analyse


    % Step 1 : find the best scoring theory for each barcode (using
    % Stouffer scores if they were calculated
  
%     % convert bestCoefs to p-values/zscores
    lengthsbar = sum(barGenRun{1}.rawBitmask);

    % intermediate step, check if this barcode is within the output struc
    % (depeends on its length)
    [barsPassThresh,locb] = ismember(cellfun(@(x) matlab.lang.makeValidName(strrep(x.name,'.tif','')),barGenRun,'un',false),bnames{m});

    if barsPassThresh~=0 && numMatchingSpecies{m}(locb)~=0
        % take theory with highest Stouffer score

        [a,bestthr] = max(arrayfun(@(x) rM{m}{x}{locb}.stoufferScore, refNums{m}{locb}));
        bT = refNums{m}{locb}(bestthr);
        pval1 = rM{m}{bT}{locb}.stoufferScore;

        
        % Step 2 : recalculate (if needed) best alignment against theory. Useful to check
        % that the alignment was correct. Not needed if rezults is saved as mat
        % file
        theoryStructSel = theoryStruct(bT); % against the same theory, but we don't know which is correct! so we calcualte it for each length separately
        % barGenRun = barGen(1);
        % w = [200:50:sum(barGenRun{1}.rawBitmask)]; % in practice could run all
        [rezOutRecalc] = local_alignment_assembly(theoryStructSel, barGenRun,mpval(m));
    
        %


        % Step 3 : block bootstrapping
        [pccScore{m} , bar1, bar2] = block_bootstrapping(barGenRun, theoryStructSel,rezOutRecalc{1},1, 2,mpval(m));
%         figure,plot(zscore(bar1,1));
%         hold on
%         plot(zscore(bar2,1))
% %     
        pccScoreV{1}{1} = cell(1,length(pccScore{m}));
        for i =1 :length(pccScore{m})
                pccScoreV{1}{1}{i} = rM{m}{bT}{locb};
                pccScoreV{1}{1}{i}.maxcoef(1) = pccScore{m}(i);
        end

        [pval] = stouffer_zscore(pccScoreV,theoryStruct,mpval(m), 0.07);
        
%         pval = arrayfun(@(x) 1-beta_ev_cdf(x, nuF*rM{m}{refNums{locb}(1)}{locb}.lengthMatch, 1, 2*max(lengthsbar, theoryStruct(refNums{locb}(1)).length),0),pccScore{m});
        match_score = cellfun(@(x) x.stoufferScore, pval{1}{1});

        % limit
        match_score(match_score>10) = 10;
        match_score(match_score<-10) = -10;

        % todo : deal with following in this case too
%         if min(pval) > 0.1 % consider only the case when pval > 0.01
%             match_score = nan;
%         else
%             pval(pval>0.9999999999999) =0.9999999999999;
%             pval(pval<10^(-20)) = 10^(-20);
%             match_score = -norminv(pval); % or norminv(1-pval)
%         end

        scores(m,:) = [mean( pccScore{m}) 3*std(pccScore{m}) length(refNums{m}{locb}) numMatchingSpecies{m}(locb) pval1 nanmean(match_score) nanstd(match_score)];
    else
    end
end

% figure,histogram(match_score)
% xlabel('match score (-norminv(pval,0,1)')
% ylabel('match score histogram')


end

