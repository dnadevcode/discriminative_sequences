%% Compare the speed of scamp vs matlab
% todo: test different platforms

pathMain = 'C:\Users\Lenovo\git\dnadevcode\'; % path to reps.

foldToRun = '';

addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))

% load some data
dirName = '/export/scratch/albertas/data_temp/bargrouping/local_results_from_tetralith/'; % work pc%
import Core.load_data_fun;
[kymoStructs,barGen,sets] = load_data_fun(dirName,5);

curDirKymos = sets.dirName; % current directory with kymographs

nmbp = sets.nmbp;

% theory loading

import Core.load_theory_structure;
thryFileIdx = 1; % todo: pass directly the theory file here
[theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);


%
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sets.theory.stretchFactors  = 0.9:0.025:1.1;


Ntheories = 10;

sets.comparisonMethod = 'mpnan';
import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen,theoryStruct(1:Ntheories), sets, [] );

%% MP: V1: concatenate theories
foldSynth=strcat('output',timestamp);
[~,~] =mkdir(foldSynth);
% have to save separately..
[namesBar, stridx,barStruct] = Core.save_bars_rescaled_txt(barGen,sF,foldSynth);
% save one long theory
[names2, baridx2] = Core.save_long_theory(theoryStruct(1:Ntheories),'barcoli');

tS{1}.filename = names2{1};

numWorkers = 30;
tic
sF = 1;
sets.w = 300;
import Core.compare_mp_multi_theories
for i=1:length(barGen)
    [compStrOut] = compare_mp_multi_theories(barGen(i),tS,sets.theory.stretchFactors,sets.w,numWorkers);
end
toc
%%

ixx=refNumsMP{wIdx}{1}(4);

bGNew = {};

bGNew{1} = barGen{1};
bGNew{2}.rawBarcode = theoryStruct(ixx).theoryBarcode;
bGNew{2}.rawBitmask = true(1,length(theoryStruct(ixx).theoryBarcode));

%     barcodeGen = bG{idxRun};
lengths = cellfun(@(x) sum(x.rawBitmask),bGNew);
tic
% bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
[oS] = calc_overlap_mp(bGNew,sF,  sets.w ,timestamp);
toc


import Core.plot_match_simple;
[f] = plot_match_simple(bGNew, oS, 2,1);


%%
foldSynth=strcat('output',timestamp);
[~,~] =mkdir(foldSynth);
% have to save separately..
[namesBar, stridx,barStruct] = Core.save_bars_rescaled_txt(barGen,sF,foldSynth);
% save one long theory
[names2, baridx2] = Core.save_long_theory(theoryStruct(1:1000),'barcoli');


%     [oS] = calc_overlap_mp(bGNew,sF,  sets.w ,timestamp);
tic
import Core.calc_overlap;
[mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,'SCAMP', sets.w, numWorkers,namesBar,names2);
toc

tic
import Core.calc_overlap_both;
[mpI1,mp1,maxMP,mpI1B,mp1B] = calc_overlap_both(barStruct,timestamp,'SCAMP', sets.w, numWorkers,namesBar,names2);

% [mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,'SCAMP', sets.w, numWorkers,names2,namesBar);
toc
% toc


for ii=1:length(baridx2)
    baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
end
% we want to create a nicer structure for all-to-all comparison and
% contains easily accesible data.
% tic
import Core.mp_res_to_struct;
[overlapStruct2] = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
% toc



%%

sets.comparisonMethod = 'mpnan';
import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen,theoryStruct(1:10), sets, [] );

numWorkers = 30;
tic
sF = 1;
sets.w = 300;
import Core.compare_mp_multi_theories
[compStrOut] = compare_mp_multi_theories(barGen,theoryStruct(1:10),sets.theory.stretchFactors,sets.w,numWorkers);
toc