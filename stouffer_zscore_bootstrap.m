function [rM] = stouffer_zscore_bootstrap(rM,  scores)
    % calculates Stouffer score

minCC = 0.5;


for i =1:length(rM)
%     i
    for j=1:length(rM{i})
       
        for k=1:length(rM{i}{j})
        % calc p-vals
        if rM{i}{j}{k}.maxcoef(1) > minCC % fix stretch factor as the smallest
            rM{i}{j}{k}.stoufferScoreNorm = (scores{k}(i,5) - rM{i}{j}{k}.stoufferScore)/scores{k}(i,6);
        else
            rM{i}{j}{k}.stoufferScoreNorm = nan;
        end
        end
    end
end



end

