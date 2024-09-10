function [cI,bI,compI2, parl,allCoefsFit,m2,mA] = get_pseudotheory_positions(bG, twoList,expPar, fastaFileF, sets,psf,sigma)

% default parameters
sF = 0.8:0.025:1.2; % length re-scaling parameter
sets.nuF = 0.08; sets.nF = 0.2; % pval params 
sets.comparisonMethod = 'mass_pcc'; % could use C++ version or UCR ED as well
sets.w = nan;
gcSF  = 1; % gc re-scaling (todo - calculate estimated number of ligands per pixel)
kY = 10; % constant yoyo
kN = 60; % constant netropsin
cN = 6; % concentration netropsin
cY = 0.02; % concentration yoyo
% pxSize = nmpx/nmbp; % pixel size
ligandLength = 4; % ligand length
isC = 1; % whether circular

if nargin < 6
psf = 370; % psf
sigma = 0.68;
end

yoyo = 26;
netrospin = 0.4;
parchangeId = 5; % scaling factor for binding constants
par = 0.68;%

mA = zeros(1,length(bG));
sA = zeros(1,size(twoList,1));
mtwo= zeros(1,size(twoList,1));
sttwo= zeros(1,size(twoList,1));
mps= zeros(1,length(bG));
allCoefsHCA = cell(1,length(bG));
allCoefsFit = cell(1,length(bG));


    hcaSets.pattern = 'binding_constant_rules.txt';
    hcaSets.method = 'literature';

bI = {};
cI = {};
hc = {};
parl = {};
import Core.run_hca_theory;

for expNr =1:length(bG)
    import CBT.Hca.Core.Theory.create_memory_struct;
    [chr1{expNr},header] = create_memory_struct(fastaFileF{expNr},num2str(expNr));
 
end

parfor expNr = 1:length(bG)
    % expNr
    barcodeGen = bG{expNr};

    % params
    nmbp = expPar{expNr}.nmbp;
    nmpx = expPar{expNr}.nmpx;
    [hcaSets,timestamp] = set_def();

    hcaSets.pattern = 'binding_constant_rules.txt';
    hcaSets.method = 'literature';

    hcaSets.pixelWidthNm = nmpx;

    pxSize = nmpx/nmbp; % pixel size
    parlist = [gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY,cN,kY,ligandLength];
    parlistcell = num2cell(parlist) ;

     hcaSets.folder{1} = chr1{expNr} ;    

    hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);
    [compI2{expNr},rezI2,m2(expNr),st2(expNr)] = run_comp(barcodeGen,hcaTheory,nmbp,sets,sF); % change this one

    
    stoufferScoresL = cellfun(@(x) double(norminv(1-x.pval)),compI2{expNr},'un',true);
    
    mAll = cellfun(@(y) cellfun(@(x) x.maxcoef(1),y), rezI2,'un',false);
    
    idx  = find(stoufferScoresL > 1.44);

    allCoefsHCA{expNr} = cellfun(@(x) x.maxcoef(1),compI2{expNr}(idx));

    mA(expNr) = mean(mAll{1}(idx));
    sA(expNr) = std(mAll{1}(idx));
    
    bI{expNr} = barcodeGen(idx);
    parl{expNr} = parlistcell;
    hc{expNr} = fastaFileF{expNr};

    [~,mid,~] = fileparts(fastaFileF{expNr});
    delete(['seq_example',mid,'_',num2str(ligandLength),'_', num2str(pxSize) ,'.mat']);

    import Helper.run_comparison_fit;
    [cI{expNr}, rezI,   mtwo(expNr) , sttwo(expNr), sS, theoryStr, yoyoBindingProb, netropsinBindingConst] = ...
        run_comparison_fit( bI{expNr},fastaFileF{expNr},parlistcell{:},parchangeId,par,sF,sets);
    allCoefsFit{expNr} = cI{expNr}.maxcoef;
end



end

