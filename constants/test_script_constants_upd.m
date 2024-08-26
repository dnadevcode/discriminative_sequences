
[hcatheory.sets,hcatheory.names] = Core.Default.read_default_sets('hcasets.txt');

hcatheorySets = hcatheory.sets.default;

hcatheorySets.folder = {'C:\Users\Lenovo\postdoc\DATA\Mapping\FASTAS\DA32087.fasta' } ;

hcatheorySets.computeBitmask = 0;
hcatheorySets.meanBpExtNm = 0.34;

% hcatheorySets%.pattern = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\new_binding_constants\netroconstants5length.txt';

disp(['N = ',num2str(length(  hcatheorySets.folder )), ' sequences to run'])

%%

yoyo = 26;
netropsin = 0.4;

ix=1
import Core.run_hca_theory;
hcaTheory = run_hca_theory(hcatheorySets,yoyo(ix),netropsin);

import CBT.Hca.UI.Helper.theory_mat_to_struct;
hcaTheory = theory_mat_to_struct(hcaTheory);


%%

nmbp = 0.21;


sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(nmbp, hcaTheory,sets );


%experiment
barcodeGen = bG{end};
sets.comparisonMethod = 'mass_pcc';
sets.w = nan;
sF = 0.9:0.025:1.1;
[compI,rezI,~] = compare_to_t(barcodeGen,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?


mean(cellfun(@(x) x.maxcoef(1),compI))


%% Run for a loop of values

barcodeGen = bG{1};
nmbp = 0.236;

yoyo =26;
netrospin = 0.4;
hcatheorySets.method = 'custom';
hcatheorySets.computeFreeConcentrations = 0;

hcatheorySets.pattern = 'binding_constant_rules.txt';
import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcatheorySets,x,netrospin),yoyo,'un',false);


import CBT.Hca.UI.Helper.theory_mat_to_struct;
hcaTheoryStruct = cellfun(@(x) theory_mat_to_struct(x),hcaTheory);


sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(nmbp, hcaTheoryStruct,sets );



%
tic
[compI,rezI,~] = compare_to_t(barcodeGen,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?
toc


m =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
st = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);





% figure,plot(yoyo,cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI))

% NEW
yoyo = 25;
netrospin = 1;
hcatheorySets.method = 'custom';
hcatheorySets.computeFreeConcentrations = 1;

hcatheorySets.pattern = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\new_binding_constants\netroconstants5.txt';
import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcatheorySets,x,netrospin),yoyo,'un',false);


import CBT.Hca.UI.Helper.theory_mat_to_struct;
hcaTheoryStruct = cellfun(@(x) theory_mat_to_struct(x),hcaTheory);


sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(nmbp, hcaTheoryStruct,sets );


tic
[compIFive,rezIFive,~] = compare_to_t(barcodeGen,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?
toc

mnew =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezIFive);
stdnew = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezIFive);


% figure,plot(yoyo,cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezIFive))


%%


figure 
title('PCC for barcodes with known ground-truth')
hold on
bar(1,[m]) 
bar(2,[mnew]) 

errorbar([m mnew],[st stdnew],'.') 
ax =gca;
ax.XTick = [1 2]
ax.XTickLabel = {'4-bg constants' , '5-bp constants'};



