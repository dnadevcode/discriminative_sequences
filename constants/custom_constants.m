
[hcaSets,timestamp] = set_def();

yoyoConst = 26;
netropsinConst = 0.4;


import CBT.Hca.Core.Theory.choose_cb_model;
[hcaSets.model ] = choose_cb_model(hcaSets.theoryGen.method,hcaSets.pattern, yoyoConst, netropsinConst);

N = 4;
 [sortedSubseq,sortv,countATs] = sorted_NT(N);
sortedDNASubseq = cellfun(@(x) int2nt(x),sortedSubseq,'UniformOutput',false)';

sortedVals = cellfun(@(x) hcaSets.model.netropsinBindingConstant(x(1),x(2),x(3),x(4)),sortedSubseq);

%% sort
sigma = 20;
constFun = arrayfun(@(x) 30*exp(-x/sigma),1:length(sortedSubseq));
figure,plot(constFun);

for i=1:length(sortedSubseq)
    hcaSets.model.netropsinBindingConstant(sortedSubseq{i}(1),sortedSubseq{i}(2),sortedSubseq{i}(3),sortedSubseq{i}(4)) = constFun(i);
end

[sortOp,sortid] =  sort(hcaSets.model.netropsinBindingConstant(:),'descend');
[ids ] = arrayfun(@(x) ind2sub([4 4 4 4],x),sortid);

% is this sorted correctly? number of ats
figure,plot(countATs(ids))
xlabel('4mers sorted by binding constant')
ylabel('AT content')
% figure,plot(sortedVals)

% [sortOp,sortid] =  sort(hcaSets.model.netropsinBindingConstant(:),'descend');
% [ids ] = arrayfun(@(x) ind2sub([4 4 4 4],x),sortid);


%% Now run analysis

barcodeGen = bG{1};
nmbp = 0.236;


% barcodeGen = bG{4};
% nmbp = 0.198;
% 
% 
% barcodeGen = bG{5};
% nmbp = 0.192;

yoyo = 20; % yoyo 26
netrospin = 0.4;
hcaSets.method = 'custom';
hcaSets.computeFreeConcentrations = 1;
hcaSets.pattern = 'binding_constant_rules.txt';

sF = 0.9:0.025:1.1;


sets.nuF = 0.08; sets.nF = 0.2; % pval params

sets.comparisonMethod = 'mass_pcc';
sets.w = nan;

sigma = [0.5];
m =[];
st = [];
compI = cell(1,length(sigma));
rezI = cell(1,length(sigma));
 hcaSets.model.netropsinBindingConstant = [];
for k=1:length(sigma)
    k

%     constFun = arrayfun(@(x) 30*exp(-x/sigma(k)),1:length(sortedSubseq));
    constFun = arrayfun(@(x) 30*exp(-x/sigma(k)),ones(1,length(sortedSubseq)).*(4-sortv));


    for i=1:length(sortedSubseq)
        hcaSets.model.netropsinBindingConstant(sortedSubseq{i}(1),sortedSubseq{i}(2),sortedSubseq{i}(3),sortedSubseq{i}(4)) = constFun(i);
    end

    tic
    import Core.run_hca_theory;
    hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin,hcaSets.model.netropsinBindingConstant),yoyo,'un',false);
    toc


%         [compI{k},rezI{k}] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF);

    [compI{k},rezI{k},m(k),st(k)] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF);
end

% stouffer scores sorted
stoufferScoresS =cellfun(@(x) double(norminv(1-x.pval)),compI{1},'un',true);

figure, histogram(cellfun(@(x) x.idx,compI{1}))

cellfun(@(x) x.idx,compI{1})
%% regular
sets.nuF = 0.08; sets.nF = 0.2; % pval params

yoyo = 26;
netrospin = 0.4;
hcaSets.pattern = 'binding_constant_rules.txt';
import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);

[compI2,rezI2,m2,st2] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF);

stoufferScoresL =cellfun(@(x) double(norminv(1-x.pval)),compI2,'un',true);

%%
sets.nuF = 0.08; sets.nF = 0.2; % pval params
yoyo = 25;
netrospin = 1;
hcaSets.method = 'custom';
hcaSets.computeFreeConcentrations = 1;
hcaSets.pattern = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\new_binding_constants\netroconstants5.txt';
import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);
[compI3,rezI3,m3,st3] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF);

% mnew =cellfun(@(y) cellfun(@(x) x.maxcoef(1),y,'un',false), rezI3,'un','false');
mnew =cellfun(@(y) cellfun(@(x) x.maxcoef(1),y,'un',true), rezI3,'un',false);

stoufferScores =cellfun(@(x) double(norminv(1-x.pval)),compI3,'un',true);

% stoufferScores =cellfun(@(x) norminv(1-x.pval),compI3,'un',true);

norminv(1-pvalLocal)

%%
ix = 4;
figure ;

title('PCC for barcodes with known ground-truth')
hold on
bar(1,[m(ix)]) 
bar(2,[m2]) 
bar(3,[m3]) 

errorbar([m(ix) m2 m3],[st(ix) st2 st3],'.') 
ax =gca;
ax.XTick = [1 2 3]
ax.XTickLabel = {'AT content-based model' , 'Literature model', '5bp'};
grid on
% ylim([0.6 0.9])
set(gca, 'YScale', 'log')


%%
figure
idxs = (stoufferScoresL>1.64);

ms1 = mean(stoufferScoresS(idxs));
ms2 = mean(stoufferScoresL(idxs));
ms3 = mean(stoufferScores(idxs));
ss1 = std(stoufferScoresS(idxs));
ss2 = std(stoufferScoresL(idxs));
ss3 = std(stoufferScores(idxs));

title('Stouffer scores for barcodes with indirect ground-truth')
hold on
bar(1,[ms1]) 
bar(2,[ms2]) 
bar(3,[ms3]) 

errorbar([ms1 ms2 ms3],[ss1 ss2 ss3],'.') 
ax =gca;
ax.XTick = [1 2 3]
ax.XTickLabel = {'AT content-based model' , 'Literature model', '5bp'};
grid on
% ylim([0.6 0.9])
set(gca, 'YScale', 'log')
% figure,plot(stoufferScoresS);hold on
% plot(stoufferScoresL)
% plot(stoufferScores)

%% 