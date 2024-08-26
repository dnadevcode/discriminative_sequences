
%% get experiment folders
dirName = '/export/scratch/albertas/download_dump/S. pyogenes all data/onlyKymo/';

import Helper.get_all_folders;
[barN, twoList] = get_all_folders(dirName);

%% get experiment closest theories
refsNames = importdata('/export/scratch/albertas/download_dump/S. pyogenes all data/ClosestTheories/refs_pyo.txt');


%% choose what exp to run
rez = cell(1,size(twoList,1));
for expNr = 1:size(twoList,1);
twoList(expNr,:)
numBars =  barN{twoList(expNr,1)}(twoList(expNr,2))

thry = refsNames{twoList(expNr,1)};

%% get theory location
allTheoryFold = '/export/scratch/albertas/download_dump/single/*.fasta';
dr = dir(allTheoryFold);

thryLoc = find(arrayfun(@(x) ~isempty(strfind(dr(x).name, thry)),1:length(dr)));

fastaFile = fullfile(dr(thryLoc).folder,dr(thryLoc).name);

%% Now load the experiments
import Helper.load_kymo_data;
[kymoStructs,barGen,nmpx,nmbp] = load_kymo_data(dirName,1,twoList(expNr,1),twoList(expNr,2));

%% Theory. All parameters

% default parameters
sets.nuF = 0.08; sets.nF = 0.2; % pval params from assembly paper
sets.comparisonMethod = 'mass_pcc'; % could use C++ version or UCR ED as well
sets.w = nan;

sF = 0.9:0.025:1.1;

sigma = 0.5;
gcSF  = 1;
kY = 10; kN = 30;
psf = 370;
cN = 6;
cY = 0.02;
pxSize = nmpx/nmbp;
isC = 1;
ligandLength = 4;
% par= 0.01:0.01:0.1;

parlist = [gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY,cN,kY,ligandLength];
parlistcell = num2cell(parlist) ;

%% Regular
[hcaSets,timestamp] = set_def();
hcaSets.folder{1} = fastaFile;

sets.nuF = 0.08; sets.nF = 0.2; % pval params

yoyo = 26;
netrospin = 0.4;
hcaSets.pattern = 'binding_constant_rules.txt';
import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);

[compI2,rezI2,m2,st2] = run_comp(barGen,hcaTheory,nmbp,sets,sF);

stoufferScoresL = cellfun(@(x) double(norminv(1-x.pval)),compI2,'un',true);

mAll = cellfun(@(y) cellfun(@(x) x.maxcoef(1),y), rezI2,'un',false);

idx  = find(stoufferScoresL > 1.44);
mA = mean(mAll{1}(idx));
sA = std(mAll{1}(idx));

%% 5mer
hcaSets.method = 'custom';
hcaSets.computeFreeConcentrations = 1;
hcaSets.pattern = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\new_binding_constants\netroconstants5.txt';
hcaSets.pattern = '/export/scratch/albertas/download_dump/netroconstants5.txt';

import Core.run_hca_theory;
hcaTheory = arrayfun(@(x) run_hca_theory(hcaSets,x,netrospin),yoyo,'un',false);
[compI3,rezI3,m3,st3] = run_comp(barGen(idx),hcaTheory,nmbp,sets,sF);

% mnew =cellfun(@(y) cellfun(@(x) x.maxcoef(1),y,'un',false), rezI3,'un','false');
mnew =cellfun(@(y) cellfun(@(x) x.maxcoef(1),y,'un',true), rezI3,'un',false);

stoufferScores =cellfun(@(x) double(norminv(1-x.pval)),compI3,'un',true);
m3
%%
% parch = {}

% run comparison, change gcSF
import Helper.run_comparison;
parchangeId = 5; % scaling factor for binding constants
par = 0.68;%
% par = [0.4:0.01:0.8];
% parchangeId = 2; % bp per px
% par = [400:10:800];
% parchangeId = 1;
% par = 0.9:0.01:1.1;
% parchangeId = 8;
% par = 0.001:0.01:0.02;
% parchangeId = 9;
% par = 1:10:200;
% parchangeId = 7; % psf
% par = 200:10:600;
[compI,rezI,m,st,sS] = run_comparison(barGen(idx),fastaFile,parlistcell,parchangeId,par,sF,sets);

% figure,plot(par,m)

rez{expNr}.m = m;
rez{expNr}.st = st;

rez{expNr}.m5 = m3;
rez{expNr}.st5 = st3;

rez{expNr}.mA = mA;
rez{expNr}.sA = sA;

end

allValues = cellfun(@(x) [x.m x.m5 x.mA],rez(1:end-2),'UniformOutput',false)
allValuesStd = cellfun(@(x) [x.st x.st5 x.sA],rez(1:end-2),'UniformOutput',false)

listOutputs = cell2mat(allValues');
%%
f=figure,bar(listOutputs)
title('Theory models comparison')
legend({'2state','Dibya (5mer)','Literature (Default in HCA)'},'location','southoutside')
xlabel('S.Pyo dataset')
ylabel('Mean centered cosine similarity')
ylim([0.5 0.8])

exportgraphics(f,'outfig.png', 'Resolution',500)


%%

figure ;

title('PCC for barcodes with known ground-truth')
hold on

bar(1,[listOutputs(:,1)]) 
bar(2,[m2]) 
bar(3,[m3]) 

errorbar([m(ix) m2 m3],[st(ix) st2 st3],'.') 
ax =gca;
ax.XTick = [1 2 3]



% 
% [a,b] = max(m);
% 
% mAllT =cellfun(@(y) cellfun(@(x) x.maxcoef(1),y), rezI{b},'un',false);
% mean(mAllT{1})
% 
% mean(mAllT{1})/mean(mAll{1}(idx))-1 % improvement in percentage of PCC score


%% Discriminative
expNr = 1;

sets.theoryFile{1} = '/export/scratch/albertas/data_temp/Alignment/data/theory/theoryGen_0.34_110_370_0_2024-08-21_16_29_08_session.mat';
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);

% 
% tic
% ticBytes(gcp);
import Helper.load_kymo_data;
[kymoStructs,barGen,nmpx,nmbp] = load_kymo_data(dirName,1,twoList(expNr,1),twoList(expNr,2));

sets.theory.nmbp = nmbp;

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct, sets );


%%

% nr = 7:8;

sets.theory.stretchFactors = sF;
sets.dirName = '/home/avesta/albertas/reps/hca/test/09_constants_test';
import Helper.calc_coefs;
[maxCoef] = calc_coefs(theoryStruct,barGen(1:5),0,sF,sets);


% unique species names
import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name({theoryStruct.name});

% discriminative locations
import Core.Discriminative.disc_locations;
[refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoef}, 0.05);

import Core.Discriminative.disc_true;
[is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);


save('maxCoef.mat','maxCoef','idx')

% uniqueSpeciesNames
uniqueSpeciesNames(idSpecies(refNums{2}))'
