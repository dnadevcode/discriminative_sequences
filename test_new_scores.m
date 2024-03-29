
% bootstrapping

scores = cell(1,length(barGen));

for ii=1:length(barGen)
    [scores{ii},pccScore] = local_bootstrap_run( barGen(ii),rM,bnames,theoryStruct ,mpval,speciesLevel,idc);

end


% todo: for all data
import Zeromodel.beta_ev_cdf;
nuF = 0.1; % should be tuned based on synthetic data? Depends a bit on pdf

import Zeromodel.beta_ev_cdf; % correct form?
newscorefun = @(x,l1,l2) norminv(beta_ev_cdf(x,nuF*l1,1,2*(max(l1,abs(l2-l1))),0));

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

f= figure,histogram(cellfun(@(x) x{8}.stoufferScore,rezMaxMP))
title('Stouffer score')
% save('')

matDirpath = regexprep(strrep(files.name,'_MP_', '_Stouffer_'), ['table_' '.*'], '');
% matDirpath = strrep(,'table_','');
import CBT.Hca.Export.export_p_vals_table;
[T] = export_p_vals_table( theoryStruct, rezMaxMP, barGen,matDirpath);


%% normalized by std

baridxs = 1:length(rezMaxMP{i});
baridxs = 1;
for i =1:length(rezMaxMP)
    i
    for j=baridxs
        % calc p-vals
        if rezMaxMP{i}{j}.maxcoef(1) > minCC 
            rezMaxMP{i}{j}.stoufferScoreNorm = (scores{j}(5) - rezMaxMP{i}{j}.stoufferScore)/scores{j}(6);
        else
            rezMaxMP{i}{j}.stoufferScoreNorm = nan;
        end
%          rezMaxMP{i}{j}.length = lengthMatch;
    end
end

figure,histogram(cellfun(@(x) x{1}.stoufferScoreNorm,rezMaxMP))
title('Stouffer score normalized by bootstrapped stDev')
