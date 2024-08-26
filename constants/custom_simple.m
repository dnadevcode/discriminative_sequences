
% load chrom data
% load('C:\Users\Lenovo\postdoc\DATA\TEMP\3all_2023-11-24_12_38_55.mat', "bG");
load('C:\Users\Lenovo\postdoc\DATA\TEMP\3_2023-11-22_11_34_44resRun.mat', "bG");
load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/3_2023-11-22_11_34_44resRun.mat', "bG");

%% Sort all 4-basepairs and count their AT ratio
N = 4;
[sortedSubseq,sortv,countATs] = sorted_NT(N);


%%
hcaSets.folder = {'C:\Users\Lenovo\postdoc\DATA\Mapping\FASTAS\DA32087.fasta' } ;
hcaSets.folder = {'/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/DA32087.fasta'};
fasta = fastaread(hcaSets.folder{1});
ntSeq = nt2int(fasta.Sequence);


% cummulative sum of AT's. 
numWsCumSum = cumsum((ntSeq == 1)  | (ntSeq == 4) );

% x = [(numWsCumSum(N+1:end) - numWsCumSum(1:end-N)) zeros(1,N-1)];% boundary condition

%%
dataIdx = 1;

barcodeGen = bG{dataIdx};

[nmbp, nmpx] = Input.extract_extension_params(fileparts(bG{dataIdx}{1}.name));

gcSF = 1; % gcsf ratio
sigma = 0.5; % at/cg conversion
netrConst = 30;
isC = 1;
pxSize = nmpx/nmbp;

[theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,netrConst,300,0.02,6,26);

 


%% compare to old
[hcaSets,timestamp] = set_def(hcaSets.folder,nmbp);

sets.nuF = 0.08; sets.nF = 0.2; % pval params from assembly paper

yoyo = 26;
netrospin = 0.4;
hcaSets.pattern = 'binding_constant_rules.txt';
import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);


%% 
figure
plot(zscore(hcaTheory{1}.theoryBarcodes{1}))
hold on
plot(zscore(theoryStr{1}.rawBarcode))

%% Compare using simple px model sigma

sF = 0.9:0.025:1.1;
sets.nuF = 0.08; sets.nF = 0.2; % pval params
sets.comparisonMethod = 'mass_pcc'; % could use C++ version or UCR ED as well
sets.w = nan;

sigma =[0.6:0.01:0.8];
gcSF = 1;
m =[];
st = [];
compI = cell(1,length(sigma));
rezI = cell(1,length(sigma));
for k=1:length(sigma)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma(k),30,300,0.02,6,20);

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    mean(stoufferScoresSimple)%(stoufferScoresSimple>1.6))
end

figure,plot(sigma,m)

%% vary gcsf

sigma = 0.68;
gcSF = 1:0.01:1.2;
sets.nuF = 0.07; sets.nF = 0.02; % pval params

m =[];
st = [];
sS = [];
compI = cell(1,length(gcSF));
rezI = cell(1,length(gcSF));
for k=1:length(gcSF)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF(k),pxSize,nmpx,isC,sigma,30,300,0.02,6,26);

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

%% Change kY
sigma = 0.68;
gcSF  = 1;
sets.nuF = 0.8; sets.nF = 0.2; % pval params
kY = [9:1:11];

m =[];
st = [];
sS = [];
compI = cell(1,length(kY));
rezI = cell(1,length(kY));
for k=1:length(kY)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,30,300,0.02,6,kY(k));

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

figure,plot(kY,m)

%% change kN
sigma = 0.68;
gcSF  = 1;
sets.nuF = 0.8; sets.nF = 0.2; % pval params
kY = 10;
kN = 5:5:50;
m =[];
st = [];
sS = [];
compI = cell(1,length(kN));
rezI = cell(1,length(kN));
for k=1:length(kN)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,kN(k),300,0.02,6,kY);

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

figure,plot(kN,m)
%% PSF
sigma = 0.68;
gcSF  = 1;
kY = 10; kN = 30;
par = 250:20:600;
m =[];
st = [];
sS = [];
compI = cell(1,length(par));
rezI = cell(1,length(par));
for k=1:length(par)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,kN,par(k),0.02,6,kY);

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

figure,plot(par,m)

%% cN

sigma = 0.68;
gcSF  = 1;
kY = 10; kN = 30;
psf = 370;
par = 1:1:10;
m =[];
st = [];
sS = [];
compI = cell(1,length(par));
rezI = cell(1,length(par));
for k=1:length(par)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,kN,psf,0.02,par(k),kY);

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

figure,plot(par,m)

%% cY

sigma = 0.68;
gcSF  = 1;
kY = 10; kN = 30;
psf = 370;
cN = 6;
par= 0.01:0.01:0.1;
m =[];
st = [];
sS = [];
compI = cell(1,length(par));
rezI = cell(1,length(par));
for k=1:length(par)
    [theoryStr] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,kN,psf,par(k),cN,kY);

    [compI,rezI,~] = compare_to_t(barcodeGen,theoryStr,sF,sets); % maybe calculate Stouffer score here also?
    m(k) =cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI);
    st(k) = cellfun(@(y) std(cellfun(@(x) x.maxcoef(1),y)), rezI);

    m(k)
    stoufferScoresSimple =cellfun(@(x) double(norminv(1-x.pval)),compI,'un',true);
    sS(k ) = mean(stoufferScoresSimple);%(stoufferScoresSimple>1.6))
end

figure,plot(par,m)