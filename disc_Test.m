
% elts = zeros(1,length(rezMaxMP)*length(rezMaxMP{1})*length(rezMaxMP{1}{1}));
% i=1;
% for idx=1:length(rezMaxMP)
%     for idy=1:length(rezMaxMP{idx})
%         for idz = 1:length(rezMaxMP{idx}{idy})
%             elts(i) = rezMaxMP{idx}{idy}{idz}.maxcoef(1)< 0.5;
%             i=i+1;
%         end
%     end
end

% cellfun(@(z) cellfun(@(y) cellfun(@(x) x{1}.maxcoef(1),y),z),rezMaxMP)
tic
    % allCoefsSingle = cellfun(@(barix) cellfun(@(y) cellfun(@(x) x.maxcoef(1),y),barix,'un',false),rezMaxMP,'un',false);
toc

% allCoefsSingle = arrayfun(@(x) arrayfun(@(y) arrayfun(@(z) rezMaxMP{z}{y}{x}.maxcoef(1) ,1:length(rezMaxMP)),1:length(rezMaxMP{1}),'un',false),1:length(rezMaxMP{1}{1}),'un',false);




import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name({theoryStruct.name});


cdiff = 0.05;
ecoliT = 1197;
theories = cell(1,size(allCoefs,1));

for barix=1:size(allCoefs,1)
    % allCoefsSingle = cellfun(@(y) cellfun(@(x) x{barix}.maxcoef(1),y),rezMaxMP,'un',false);
    
    curBdisc = zeros(1,size(allCoefs,2));
    for barid = 1:size(allCoefs,2)
        singleCoefs =  squeeze(allCoefs(barix,barid,:));
        
        [a,sortedid] = sort(singleCoefs,'desc','MissingPlacement','last');
        %         maxs = a(1);
        
        discLocations = (find(a>a(1)-cdiff));
        theories{barix}{barid} = sortedid(discLocations);
    end
end
    
is_distinct = [];
uniqueMatchSequences = [];
for i=1:length(theories)
    import Core.Discriminative.disc_true;
    [is_distinct{i}, numMatchingSpecies, uniqueMatchSequences{i}] = disc_true(theories{i}, idSpecies);
end


figure,imagesc(cell2mat(is_distinct'));colormap gray
xlabel('Barcode nr.')
ylabel('Re-scaling factor nr.')

%     idSpecies(cellfun(@(x) x(1),refNums{i})==1197)

        numPyoFirst(i) = sum(idSpecies(cellfun(@(x) x(1),refNums{i}))==1197);

    curBdisc(barid) = sum(idSpecies(theories)~=ecoliT) == 0;
% end
%%

%curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/pyo_local/S.pyogenes 388 2023-07-13 kymos/';

% curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/test/test3/km';
% curDirKymos = '/proj/snic2022-5-384/users/x_albdv/data/pyo_local/dist/kymos/';
% import Core.load_local_alignment_results_from_files;
% [rM, bnames, mpval,thryNames,files] = load_local_alignment_results_from_files(curDirKymos ); 

% probably the same for all?

% 
% sF = zeros(1,length(files));
% for i=1:length(files)
%     setsSaved = load(fullfile(files(i).folder,[strrep(files(i).name(1:end-4),'sample1_table_','run_settings'),'.mat']));
%     sF(i) = setsSaved.sets.sF;
% end
% 
%         for j=1:length(rM{2})
%             for k=1:length(rM{2}{1})
%                 rM{2}{j}{k}.maxcoef = rM{2}{j}{k}.maxcoef/5;
%             end
% 
%         end
% 
% is_distinct = cell(1,length(rM));
% numPyoFirst = zeros(1,length(rM));
% refNums =  cell(1,length(rM));
% uniqueMatchSequences = cell(1,length(rM));
% for i=1:length(rM)
%     % if i==1
%     % 
%     % end
% 
% 
%     import Core.disc_locs;
%     [refNums{i}, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{i},0,0.8);
% 
%         [refNums{i}, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{i},0,0.05);
% 
%     import Core.Discriminative.disc_true;
%     [is_distinct{i},numMatchingSpecies,uniqueMatchSequences{i}] = disc_true(refNums{i}, idSpecies);
% 
%     % uniqueMatchSequences == 1197
%     % idSpecies(cellfun(@(x) x(1),refNums{i})==1197)
%     % number of sequences with pyo as first
%     % numPyoFirst(i) = sum(idSpecies(cellfun(@(x) x(1),refNums{i}))==1197);
% end
% 
% 
% numD = cellfun(@(x) sum(x), is_distinct);
% 
% [svals,sidx] = sort(sF);
% 
% %% Figure for meeting
% figure,tiledlayout(2,1)
% nexttile
% plot(svals,numD(sidx),'o-')
% xlabel('Length re-scaling factor %')
% ylabel('Num discriminative barcodes')
% 
% nexttile
% plot(svals,numPyoFirst(sidx),'o-')
% xlabel('Length re-scaling factor %')
% ylabel('Num S.Pyogenes best match')
% 
% %% We can show the best match for one of the examples which is not significant
% distEl1 = find(is_distinct{4}==1);
% distEl2 = find(is_distinct{1}==1);
% becomeNonDistinctOneTwo = find(ismember(distEl1,distEl2)==0);
% becomeNonDistinctTwoOne = find(ismember(distEl2,distEl1)==0);
% 
% barid = distEl2(becomeNonDistinctTwoOne);
% 
% rM{3}{refNums{3}{barid(1)}(1)}{barid(1)}.or = 1; % this needs to be save in table..
% super_quick_plot_rezmax(barid(1),barcodeGenC,rM{3},theoryStruct,refNums{3}{barid(1)}(1))
% 
% thrIdxs = idSpecies(refNums{1}{barid(1)});
% badThr = find(thrIdxs~=3284,1,'first');
% ixT = refNums{i}{barid(1)}(badThr);
% 
% rM{1}{ixT}{barid(1)}.or = -1; % this needs to be save in table..
% super_quick_plot_rezmax(barid(1),barcodeGenC,rM{1},theoryStruct,ixT)
% 
% rM{3}{ixT}{barid(1)}.or = -1; % this needs to be save in table..
% super_quick_plot_rezmax(barid(1),barcodeGenC,rM{3},theoryStruct,ixT)
% 
% 
% %
%     [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
%         calculate_sorted_pvals(oS,sets.minOverlap, sets.nupar,sets.pthresh);
% %     [sortedVals, sortedIds,pscores,fullScores,overlaplen,partialScore,partialLength,pval,pvalLeftOver,pvalCombined] = sorted_scores(oS);
% import Core.barcode_island_output;
% goodPos = sortedVals > 3;
% sortedIdsGood = sortedIds(find(goodPos));
% sortedValsGood = sortedVals(find(goodPos));
% [barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGen,barcodeIslandsData, barStruct,barIslands]=...
%     barcode_island_output(sortedValsGood,sortedIdsGood, oS, bars,'tt',sets.scDiffSetting,sets.pxDifSetting, [],1);
% 
% 
% idx = [2 3]%barIslands{5};
% loc = zeros(1,length(idx));
% for i=1:length(idx)
%     [~,nm,~] = fileparts(bars{idx(i)}.name);
%     
%     nameToFind = matlab.lang.makeValidName(nm);
%     loc(i) = find(cellfun(@(x) ~isempty(strfind(x, nameToFind)),bnames{1}));
% 
% end
%  
% is_distinct(loc)
% 
% numMatchingSpecies(loc)
% 
% % What if we keep only distinct ones
