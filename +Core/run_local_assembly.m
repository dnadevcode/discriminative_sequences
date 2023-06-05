function [outputArg1,outputArg2] = run_local_assembly(barGen,sF, minOverlap,timestamp)
sF = 0.95:0.01:1.05;
minOverlap = 200;
% import Core.calc_overlap_mp
[oS] = calc_overlap_mp(barGen,sF, minOverlap,timestamp);

barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
% 
idxPair =   3;
[sortedVals, sortedIds,pscores,fullScores] = sorted_scores(oS);
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, oS,curSink,curSource);

[f] = plot_match_simple(barStruct, oS,16,19);

%%

minLen = [200:50:2000];
[MP,mpMax,theoryStructRev,MPI] = ...
bargrouping_minimum_length([],130,0.25,1,30,minLen, 5*10^6); % if thry not
% know
% [MP,mpMax,theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,110,0.25,1,30,minLen, 0);

thresCC = 0.91;
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

end

