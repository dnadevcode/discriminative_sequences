function [cI,parI, fastaN, goodScores, allScores] = get_pseudotheory_best_scores_two_state(bG, expPar, fastaFileF, sets,psf,sigma, sF)
% this provides best scores but only for the literature method

% default parameters
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

% keep this unchanged
parchangeId = 5; % scaling factor for binding constants
par = sigma;

% chr1 = cell(1,length(bG));
% for expNr =1:length(bG)
%     [chr1{expNr},~] = create_memory_struct(fastaFileF{expNr},num2str(expNr));
% end

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

% alphaNu = sets.nuF;
% alphaN = sets.nF;
% pvalfun = @(x,l1,l2) 1-beta_ev_cdf(x,alphaNu*l1,1,alphaN*2*l2,1);

import Helper.run_comparison_fit;


parfor expNr = 1:length(bG)
    % expNr
    barcodeGen = bG{expNr};
    % params take from experiment
    nmbp = expPar{expNr}.nmbp;
    nmpx = expPar{expNr}.nmpx;
    [hcaSets, timestamp] = set_def();
    hcaSets.comparisonMethod = 'mass_pcc';
    hcaSets.w = 0;
    hcaSets.pixelWidthNm = nmpx;
    hcaSets.theory.stretchFactors = sF;

    pxSize = nmpx/nmbp; % pixel size
    parlist = [gcSF, pxSize, nmpx, isC, sigma, kN, psf, cY, cN, kY, ligandLength];
    parlistcell = num2cell(parlist) ;

    [~,mid,~] = fileparts(fastaFileF{expNr});
    delete(['seq_example',mid,'_',num2str(ligandLength),'_', num2str(pxSize) ,'.mat']);

    [cI{expNr}] = ...
        run_comparison_fit( barcodeGen,fastaFileF{expNr}, parlistcell{:}, parchangeId, par, sF, sets);

    goodScores(expNr) = mean(cI{expNr}.maxcoef);
    allScores{expNr} = cI{expNr}.maxcoef
%     cI{expNr} = compI;
    parI{expNr} = parlistcell;
    fastaN{expNr} = fastaFileF{expNr};

    %%
end



end

