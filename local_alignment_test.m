% prep_thry_for_local_comp calculates theories.
%% Need to generate proper thry
load('theories.mat');
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%%

sets.keep=1;

files = dir('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/Ecoli_DA65808/RawKymographs 1908XX - 110nm - 0.267nm_bp/*.tif');
sets.kymosets.filenames = arrayfun(@(x) files(x).name,1:length(files),'un',false);
sets.kymosets.kymofilefold = arrayfun(@(x) files(x).folder,1:length(files),'un',false);

simdata = 0;
sets.output.matDirpath = 'output';
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 6;
sets.minLen = 300;
sets.genConsensus  = 0;

% import CBT.Hca.Import.add_kymographs_fun;
% [kymoStructs] = add_kymographs_fun(sets);
   


% we have theoryStruct. Now want to compare single barcode



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
%% ASSEMBLY
sF = 0.95:0.01:1.05;
minOverlap = 200;
[oS] = calc_overlap_mp(barGen,sF, minOverlap,timestamp);





%%

tic
sets.theory.stretchFactors = 1;
sets.w = 300;
sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
% compare theory to experiment
import CBT.Hca.Core.Comparison.compare_distance;
[rezMax,bestBarStretch,bestLength] = compare_distance(barGen,theoryStruct, sets, [] );
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
[rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen,theoryStruct, sets, [] );
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

%%
barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
% 
idxPair =   20;
[sortedVals, sortedIds,pscores,fullScores] = sorted_scores(oS);
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, oS,curSink,curSource);

%%

minLen = [200:50:2000];
[MP,mpMax,theoryStructRev,MPI] = ...
bargrouping_minimum_length([],110,0.25,1,30,minLen, 5*10^6); % if thry not
% know
% [MP,mpMax,theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,110,0.25,1,30,minLen, 0);

thresCC = 0.8;
import Core.filter_good_scores
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMax,minLen,thresCC);

% plot the graph
import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIds,length(pscores),barGen,zeros(1,length(pscores)),100,oS);

%% 
import Core.create_barset_os; 
[barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = ...
    create_barset_os(sortedVals,sortedIds, barStruct,oS);

nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands)); % keep only non-empty
barcodeIslandsData = barcodeIslandsData(nonEmpty);
barcodeIslands = barcodeIslands(nonEmpty);
barislandLength = cellfun(@(x) max(x(:,2))-min(x(:,1)),barcodeIslandsData)


% plot group quick
% sourceGroup = 2;
% plot_temp(barcodeIslandsData{sourceGroup}, barcodeIslands{sourceGroup},barStruct )


% assembly w/o reference / bars are shuffled..
nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
import Plot.plot_best_pos;
fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6


[outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(barGen(cellfun(@(x) x.barid,posStruct))',posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp,1);


import Plot.plot_island;
plot_island(posStruct,islandElts, barGen(cellfun(@(x) x.barid,posStruct))',1);
