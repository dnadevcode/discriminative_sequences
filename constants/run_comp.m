function [compI,rezI,m,st] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF)

import CBT.Hca.UI.Helper.theory_mat_to_struct;
hcaTheoryStruct = cellfun(@(x) theory_mat_to_struct(x),hcaTheory);


sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(nmbp, hcaTheoryStruct,sets );



%
tic
[compI,rezI,~] = compare_to_t(barcodeGen,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?
toc


try
m =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
st = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);
catch
m =0;
st = 0;
end

end

