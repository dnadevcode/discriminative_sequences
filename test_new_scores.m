% todo: for all data
import Zeromodel.beta_ev_cdf;
nuF = 0.1; % should be tuned based on synthetic data? Depends a bit on pdf

import Zeromodel.beta_ev_cdf; % correct form?
newscorefun = @(x,l1,l2) norminv(beta_ev_cdf(x,nuF*l1,1,2*(max(l1,l2-l1)),0));

minCC = 0.5;


for i =1:length(rezMaxMP)
    i
    for j=1:length(rezMaxMP{i})
        % calc p-vals
        if rezMaxMP{i}{j}.maxcoef(1) > minCC 
            rezMaxMP{i}{j}.stoufferScore = newscorefun(rezMaxMP{i}{j}.maxcoef(1),rezMaxMP{i}{j}.lengthMatch,theoryStruct(i).length);
        else
            rezMaxMP{i}{j}.stoufferScore = nan;
        end
%          rezMaxMP{i}{j}.length = lengthMatch;
    end
end

% figure,histogram(cellfun(@(x) x{8}.stoufferScore,rezMaxMP))

matDirpath = regexprep(strrep(files.name,'_MP_', '_Stouffer_'), ['table_' '.*'], '');
% matDirpath = strrep(,'table_','');
import CBT.Hca.Export.export_p_vals_table;
[T] = export_p_vals_table( theoryStruct, rezMaxMP, barGen,matDirpath);
