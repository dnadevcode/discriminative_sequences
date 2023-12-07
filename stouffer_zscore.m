function [rM] = stouffer_zscore(rM,theoryStruct,mpval, nuF)
    % calculates Stouffer score


% todo: for all data
import Zeromodel.beta_ev_cdf;
if nargin < 4
    nuF = 0.1; % should be tuned based on synthetic data? Depends a bit on pdf
end

import Zeromodel.beta_ev_cdf; % correct form?
% newscorefun = @(x,l1,l2,w) norminv(beta_ev_cdf(x,nuF*w,1,max(0.2*2*(max(l1,l2)-w),0.01*2*(l1-w+1)*(l2-w+1)),0));
newscorefun = @(x,l1,l2,w) norminv(beta_ev_cdf(x,nuF*w,1,max(0.2*2*(max(l1,l2)-w),0.01*2*(l1-w+1)*(l2-w+1)),0));

newscorefun2 = @(x,l1,l2,w) norminv(beta_ev_cdf(x,nuF*l1,1,0.2*2*l2,0));


minCC = 0.5;


parfor i =1:length(rM)
%     i
    if mpval(i)==0
        fun = newscorefun2;
    else
        fun = newscorefun;
    end

    for j=1:length(rM{i})
%         j
        for k=1:length(rM{i}{j})
        % calc p-vals
        if rM{i}{j}{k}.maxcoef(1) > minCC % fix stretch factor as the smallest
            rM{i}{j}{k}.stoufferScore = fun(rM{i}{j}{k}.maxcoef(1),rM{i}{j}{k}.lengthMatch/rM{i}{j}{k}.bestBarStretch*0.9, theoryStruct(i).length,rM{i}{j}{k}.lengthMatch);
%             rM{i}{j}{k}.pval = fun(rM{i}{j}{k}.maxcoef(1),rM{1}{j}{k}.lengthMatch/rM{1}{j}{k}.bestBarStretch*0.9, theoryStruct(i).length,rM{i}{j}{k}.lengthMatch);

        else
            rM{i}{j}{k}.stoufferScore = nan;
        end
%          rezMaxMP{i}{j}.length = lengthMatch;
        end
    end
end



end

