function [] = run_analysis_5()


N = 5;
 [sortedSubseq,sortv,countATs] = sorted_NT(N);


barcodeGen = bG{1};
nmbp = 0.236; %todo: extract from name

yoyo = 20; % yoyo 26
netrospin = 0.4;
hcaSets.method = 'custom';
hcaSets.computeFreeConcentrations = 0;
hcaSets.pattern = 'binding_constant_rules.txt';

sF = 0.9:0.025:1.1;


sets.nuF = 0.08; sets.nF = 0.2; % pval params

sets.comparisonMethod = 'mass_pcc';
sets.w = nan;

sigma = 1;
m =[];
st = [];
compI = cell(1,length(sigma));
rezI = cell(1,length(sigma));
hcaSets.model.netropsinBindingConstant = zeros(4,4,4,4,4);
for k=1:length(sigma)
    k

%     constFun = arrayfun(@(x) 30*exp(-x/sigma(k)),1:length(sortedSubseq));
    constFun = arrayfun(@(x) 30*exp(-x/sigma(k)),ones(1,length(sortedSubseq)).*(5-sortv));


    for i=1:length(sortedSubseq)
        hcaSets.model.netropsinBindingConstant(sortedSubseq{i}(1),sortedSubseq{i}(2),sortedSubseq{i}(3),sortedSubseq{i}(4),sortedSubseq{i}(5)) = constFun(i);
    end

    tic
    import Core.run_hca_theory;
    hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin,hcaSets.model.netropsinBindingConstant),yoyo,'un',false);
    toc

     tic
    import Core.run_hca_theory;
    hcaTheory2 = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);
    toc

    [compI{k},rezI{k},m(k),st(k)] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF);
end

% stouffer scores sorted
stoufferScoresS =cellfun(@(x) double(norminv(1-x.pval)),compI{1},'un',true);

figure, histogram(cellfun(@(x) x.idx,compI{1}))

cellfun(@(x) x.idx,compI{1})


end

