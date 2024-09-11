function [bI, cI,parI, fastaN, goodScores, allScores] = get_pseudotheory_best_scores_literature(bG, twoList, expPar, fastaFileF, sets,psf,sigma, sF)
% this provides best scores but only for the literature method

% default parameters // most are unused
sets.nuF = 0.08; sets.nF = 0.2; % pval params
sets.comparisonMethod = 'mass_pcc'; % could use C++ version or UCR ED as well
sets.w = nan;
gcSF  = 1; % gc re-scaling (todo - calculate estimated number of ligands per pixel)
kY = 10; % constant yoyo
kN = 60; % constant netropsin
cN = 6; % concentration netropsin
cY = 0.02; % concentration yoyo
ligandLength = 4; % ligand length
isC = 1; % whether circular
yoyo = 26;
netrospin = 0.4;

if nargin < 6
    psf = 370; % psf
    sigma = 0.68;
    sF = 0.8:0.025:1.2; % length re-scaling parameter
end

chr1 = cell(1,length(bG));
for expNr =1:length(bG)
    [chr1{expNr},~] = create_memory_struct(fastaFileF{expNr},num2str(expNr));
end

goodScores = zeros(1,length(bG));
% sA = zeros(1,size(twoList,1));
% allCoefsHCA = cell(1,length(bG));
% 
% hcaSets.pattern = 'binding_constant_rules.txt';
% hcaSets.method = 'literature';

bI = cell(1,length(bG)); % only good barcodes
cI = cell(1,length(bG));
fastaN = {};
parStruct = {};
import Core.run_hca_theory;

import CBT.Hca.Core.Theory.create_memory_struct;
import CBT.Hca.Core.Comparison.hca_compare_distance;
import Core.Discriminative.generate_sf_struct;

import Zeromodel.beta_ev_cdf; % pcc

alphaNu = sets.nuF;
alphaN = sets.nF;
pvalfun = @(x,l1,l2) 1-beta_ev_cdf(x,alphaNu*l1,1,alphaN*2*l2,1);


parfor expNr = 1:length(bG)
    % expNr
    barcodeGen = bG{expNr};
    % params take from experiment
    nmbp = expPar{expNr}.nmbp;
    nmpx = expPar{expNr}.nmpx;
    [hcaSets, timestamp] = set_def();

    hcaSets.pattern = 'binding_constant_rules.txt';
    hcaSets.method = 'literature';
    hcaSets.comparisonMethod = 'mass_pcc';
    hcaSets.w = 0;
    hcaSets.pixelWidthNm = nmpx;
    hcaSets.theory.stretchFactors = sF;

    pxSize = nmpx/nmbp; % pixel size
    parlist = [gcSF, pxSize, nmpx, isC, sigma, kN, psf, cY, cN, kY, ligandLength];
    parlistcell = num2cell(parlist) ;

    hcaSets.folder{1} = chr1{expNr} ;

    hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets, x, netrospin), yoyo,'un',false);

    import CBT.Hca.UI.Helper.theory_mat_to_struct;
    hcaTheoryStruct = cellfun(@(x) theory_mat_to_struct(x),hcaTheory);

    hcaSets.theory.theoryDontSaveTxts = 1;
    import CBT.Hca.Core.Analysis.convert_nm_ratio;
    theoryStruct = convert_nm_ratio(nmbp, hcaTheoryStruct,hcaSets );

    [rezMax] = hca_compare_distance(barcodeGen, theoryStruct, hcaSets );

    [compI] = generate_sf_struct(rezMax,hcaSets);

    % only for masspcc
    compI.pval = pvalfun(compI.maxcoef,double(compI.bestlength),theoryStruct.length);

    compI.stoufferScoresSimple = double(norminv(1 - compI.pval));
%     sS(k ) = mean( compI.stoufferScoresSimple); %(stoufferScoresSimple>1.6))

    idx  = find(compI.stoufferScoresSimple  > 1.44);

    bI{expNr} = barcodeGen(idx);
    goodScores(expNr) = mean(compI.maxcoef(idx));
    allScores{expNr} = compI.maxcoef(idx);
    cI{expNr} = compI;
    parI{expNr} = parlistcell;
    fastaN{expNr} = fastaFileF{expNr};

    %%

%     [scoreStruct{expNr},rezI2,m2(expNr),st2(expNr)] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF); % change this one


%     stoufferScoresL = cellfun(@(x) double(norminv(1-x.pval)),scoreStruct{expNr},'un',true);

%     allScores = cellfun(@(y) cellfun(@(x) x.maxcoef(1),y), rezI2,'un',false);

%     idx  = find(stoufferScoresL > 1.44);
% 
%     allCoefsHCA{expNr} = cellfun(@(x) x.maxcoef(1),scoreStruct{expNr}(idx));
% 
%     goodScores(expNr) = mean(allScores{1}(idx));
%     sA(expNr) = std(allScores{1}(idx));

%     parStruct{expNr} = parlistcell;
%     hc{expNr} = fastaFileF{expNr};
end



end

