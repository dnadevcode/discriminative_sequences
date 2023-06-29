% prep_thry_for_local_comp calculates theories.
load('theories.mat');

% we have theoryStruct. Now want to compare single barcode

sets.keep=1;
sets.kymosets.filenames ={ '20221220_87-st7_filter-2_int-45_mol-55_molecule_1_kymograph.tif'};
sets.kymosets.kymofilefold = {'/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/ecoli_2_new/12_20_2/_sessiondata/'};
simdata = 0;
sets.output.matDirpath = 'output';
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 6;
sets.minLen = 300;
sets.genConsensus  = 0;

if simdata
    tic
    N = 500;
    A = imgaussfilt(normrnd(0,1,1,N),3);
    
    M = 50;
    barGen = [];
    for i=1:M
        barGen{i}.rawBarcode = A;
        barGen{i}.rawBitmask = ones(1,length(A));
    end
else
    import CBT.Hca.Import.add_kymographs_fun;
    [kymoStructs] = add_kymographs_fun(sets);
    
    sets.alignMethod = 1;
    sets.edgeDetectionSettings.method = 'Zscore';
    % align kymos
    import CBT.Hca.Core.align_kymos;
    [kymoStructs] = align_kymos(sets,kymoStructs);
    
    % generate barcodes
    import CBT.Hca.Core.gen_barcodes;
    barGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);


end
%%

tic
sets.theory.stretchFactors = 1;
sets.w = 300;
sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
% compare theory to experiment
import CBT.Hca.Core.Comparison.compare_distance;
[rezMax,bestBarStretch,bestLength] = compare_distance(barGen,theoryStructT, sets, [] );
%
barIx = 2;
toc
allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMax);

% figure,plot(allCCs)
% xlabel('Theory bar')
% ylabel('Score')

[a,sortedid] = sort(allCCs,'desc','MissingPlacement','last');
maxs = a(1);
cdiff = 0.05;
selectedRef = sortedid(find(a>a(1)-cdiff))

figure,plot(a(1:10))
hold on
plot([0 10],[maxs-cdiff maxs-cdiff])
xlabel('Theory bar')
ylabel('Score')
legend({'Score','Cut-off'})

%% mpnan
tic
sets.w = 400;
sets.comparisonMethod = 'mpnan';

import CBT.Hca.Core.Comparison.compare_distance;
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen,theoryStructT, sets, [] );
%%%


toc
allCCsMP = cellfun(@(x) x{1}.maxcoef(1),rezMaxMP);

% figure,plot(allCCsMP)
% xlabel('Theory bar')
% ylabel('Score MP')
%
[aMP,sortedidMP] = sort(allCCsMP,'desc','MissingPlacement','last');
maxsMP = aMP(1);
cdiff = 0.05;
selectedRefMP = sortedidMP(find(aMP>aMP(1)-cdiff))

theoryStruct{selectedRefMP(1)}

figure,plot(aMP(1:100))
hold on
plot([0 100],[maxsMP-cdiff maxsMP-cdiff])
xlabel('Theory bar')
ylabel('Score')
legend({'Score','Cut-off'})


%%
name = strsplit(theoryStruct{selectedRef(2)}.filename,'theory_');
nameSplit = strsplit(name{2},'_');
seqName = [nameSplit{1} '_' nameSplit{2}];

% seqName = 'NZ_CP019000';
find(cellfun(@(x) contains(x,seqName),plasmidsDataTable.RefSeq))


[plasmidsDataTable, originalColHeaders] = NCBI.get_plasmids_report();
isBacterialPlasmid = ismember(plasmidsDataTable{:,'Kingdom'}, {'Bacteria'});
dataTableBacterialPlasmids = plasmidsDataTable(isBacterialPlasmid,:);
 
%   find(cellfun(@(x) ~isempty(strfind(seqName,x)),plasmidsDataTable.RefSeq))

%% local matchings scamp
numWorkers = 30;
tic
sF = 1;
import Core.compare_mp_multi_theories
[compStrOut] = compare_mp_multi_theories(barGen,theoryStruct(1:10),sets.theory.stretchFactors,sets.w,numWorkers);
toc

% [compStr] = compare_to_t_mp(barGen,theoryStruct,sets.theory.stretchFactors,sets.w)
% toc

%% local to global


% comparisonStructAll = rezMax;
% for i=1:length(comparisonStructAll)
%     for j=1:length(bestBarStretch{i})
%         comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
%         comparisonStructAll{i}{j}.length = bestLength{i}(j);
%     end
% end
% import CBT.Hca.Core.Comparison.combine_theory_results;
% [comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);
% 
% maxcoef = comparisonStruct{1}.maxcoef;
% f=figure;
% tiledlayout(1,1)
% ax = nexttile
% import CBT.Hca.UI.Helper.plot_best_bar_mp;
% resultStruct=plot_best_bar_mp(ax,barGen,[],comparisonStruct, theoryStruct, maxcoef,1,sets);
% %

%%
% todo: expand local to full
rezMaxMP{1}{1}
thr = importdata(theoryStruct{1}.filename);
theoryStruct{1}

% recover score
bar2= thr(rezMaxMP{1}{1}.secondPos(1):rezMaxMP{1}{1}.secondPos(1)+rezMaxMP{1}{1}.lengthMatch-1);
bar1 = barGen{1}.rawBarcode(rezMaxMP{1}{1}.secondPos(1)-rezMaxMP{1}{1}.pos(1)+1:rezMaxMP{1}{1}.secondPos(1)-rezMaxMP{1}{1}.pos(1)+rezMaxMP{1}{1}.lengthMatch);

pcc = zscore(bar1,1)*zscore(bar2',1)/length(bar1)

%% full score:
posdif = rezMaxMP{1}{1}.secondPos(1)-rezMaxMP{1}{1}.pos(1)+1;
w = rezMaxMP{1}{1}.lengthMatch;
bar1 = barGen{1}.rawBarcode;
bar2 =  thr(rezMaxMP{1}{1}.secondPos-posdif+1:rezMaxMP{1}{1}.secondPos-posdif+length(bar1));

bar2=bar2(barGen{1}.rawBitmask)
bar1 = bar1(barGen{1}.rawBitmask);
pcc = zscore(bar1,1)*zscore(bar2',1)/length(bar1)

bar1 = barGen{1}.rawBarcode;
bar1out = bar1(barGen{1}.rawBitmask);

barcodeGen
full_score = zeros(1,length(rezMaxMP));
barcodeGenThry = barcodeGen(cellfun(@(x) ~isempty(x.rawBarcode),barcodeGen))
for i=1:length(rezMaxMP)
    i
%     thr = importdata(theoryStruct{i}.filename);
    thr = barcodeGenThry{i}.rawBarcode;
    thr2 = [thr thr];
    posdif = rezMaxMP{i}{1}.secondPos(1)-rezMaxMP{i}{1}.pos(1)+1;
    w = rezMaxMP{i}{1}.lengthMatch;
    try
    bar2 =  thr2(rezMaxMP{i}{1}.secondPos-posdif+1:rezMaxMP{i}{1}.secondPos-posdif+length(bar1));
    catch
        bar2 =  thr2(rezMaxMP{i}{1}.secondPos-posdif+1+length(thr):rezMaxMP{i}{1}.secondPos-posdif+length(bar1)+length(thr));
    end
    bar2out=bar2(barGen{1}.rawBitmask);
    full_score(i) = zscore(bar1,1)*zscore(bar2',1)/length(bar1);

end

%%

[aOverlap,sortedFull] = sort(full_score,'desc','MissingPlacement','last');
maxsMP = aOverlap(1);
cdiff = 0.05;
selectedfull = sortedFull(find(aOverlap>aOverlap(1)-cdiff))

theoryStruct{selectedfull(1)}

figure,plot(aOverlap(1:10))
hold on
plot([0 10],[maxsMP-cdiff maxsMP-cdiff])
xlabel('Theory bar')
ylabel('Score')
legend({'Score','Cut-off'})


%% Local - how many significant matches
refNums = cell(1,length(barGen));
for barIx = 1:length(barGen);
%     toc
    allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMax);

% figure,plot(allCCs)
% xlabel('Theory bar')
% ylabel('Score')

    [a,sortedid] = sort(allCCs,'desc','MissingPlacement','last');
    maxs = a(1);
    cdiff = 0.05;
    selectedRef = sortedid(find(a>a(1)-cdiff));
    refNums{barIx} = selectedRef;
end
