function [scoresStouffer] = stouffer_zscore_mat(scores, lengthsBars, lengthsThry, mpval, nuF,wOverlap)
    % calculates Stouffer score
    % todo: also for mp (add w)


    if nargin < 6
        wOverlap = nan;
    end

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

scoresStouffer = nan(size(scores));

if mpval==0
    fun = newscorefun2;
else
    fun = newscorefun;
end
for i =1:size(scores,1)
    %     i


    for j=1:size(scores,2)
        %         j
        lenE = lengthsBars(j);
        lenT =  lengthsThry(j);
        parfor k=1:size(scores,3)
            % calc p-vals
            if scores(i,j,k) > minCC % fix stretch factor as the smallest
                scoresStouffer(i,j,k) = fun(scores(i,j,k),lenE ,lenT, wOverlap);
            end
        end
    end
end



end

