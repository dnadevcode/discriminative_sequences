foldE = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/';
% foldE = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local_pyo/';

cont = dir(foldE)
cont = cont(and([cont.isdir], ~ismember({cont.name}, {'.', '..'})));

coefs = cell(1,length(cont));
scoresStouffer = cell(1,length(cont));
for i = 1:length(cont)
    file = dir(fullfile(fullfile(cont(i).folder,cont(i).name),'*sf_rez*'));
    coefs{i} = load(fullfile(file(1).folder,file(1).name));

    % thse could be saved separately as mat files
    filebars =  dir(fullfile(fullfile(cont(i).folder,cont(i).name),'*bars*'));
    bars = load(fullfile(filebars(1).folder,filebars(1).name));
    % should save these as well..
    [scoresStouffer{i}.allCoefs] = Core.stouffer_zscore_mat(coefs{i}.allCoefs, cellfun(@(x) sum(x.rawBitmask),bars.barGen), [theoryStruct.length].*theoryStruct(1).meanBpExt_nm/0.2, 0, 0.07);

end

% we can calculate Stouffer coefficients for bars & theoryStruct
% [rM] = stouffer_zscore(rM,theoryStruct,0, 0.07)

import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name({theoryStruct.name});


%% Cascading discriminative barcodes
cdiff = 0.05;
sf  = coefs{1}.sets.theory.stretchFactors;
numCoefs =  floor(length(sf)/2)+1;

sfCascade = arrayfun(@(x) numCoefs-x:numCoefs+x,0:numCoefs-1,'un',false);

%% cascading
theories = cell(1,length(coefs));%,length(sfCascade),size(coefs{1}.allCoefs,2));
for i=1:length(coefs) % for each dataset (later move to a function)
    for j=1:length(sfCascade)
        for barid = 1:size(coefs{i}.allCoefs,2)
            singleCoefs =  squeeze(max(coefs{i}.allCoefs(sfCascade{j},barid,:),[],1));
            [a,sortedid] = sort(singleCoefs,'desc','MissingPlacement','last');
        
            discLocations = (find(a>a(1)-cdiff));
            theories{i}{j}{barid} = sortedid(discLocations);
        end
    end
end

corSp = 1197;
% corSp = 3284;

corSpVecCasc = zeros(1,length(sfCascade));
is_distinct_casc = cell(1,length(coefs));
uniqueMatchSequences = cell(1,length(coefs));
for i=1:length(coefs) % for each dataset (later move to a function)
    is_distinct_casc{i} = zeros(length(sfCascade),size(coefs{i}.allCoefs,2));
    for j=1:length(sfCascade)
        import Core.Discriminative.disc_true;
        [is_distinct_casc{i}(j,:), numMatchingSpecies, uniqueMatchSequences{i}{j}] = disc_true(theories{i}{j}, idSpecies);
        sum(cell2mat(uniqueMatchSequences{i}{j}(find(is_distinct_casc{i}(j,:))))~=corSp)
        corSpVecCasc(j) = corSpVecCasc(j) + sum(cell2mat(uniqueMatchSequences{i}{j}(find(is_distinct_casc{i}(j,:))))~=corSp);

    end
end

%% non-cascading

theoriesNoncasc = cell(1,length(coefs));%,length(sfCascade),size(coefs{1}.allCoefs,2));
for i=1:length(coefs) % for each dataset (later move to a function)
    for j=1:length(sf)
        for barid = 1:size(coefs{i}.allCoefs,2)
            singleCoefs =  squeeze(coefs{i}.allCoefs(j,barid,:));
            [a,sortedid] = sort(singleCoefs,'desc','MissingPlacement','last');
        
            discLocations = (find(a>a(1)-cdiff));
            theoriesNoncasc{i}{j}{barid} = sortedid(discLocations);
        end
    end
end



is_distinct_noncasc = cell(1,length(coefs));
uniqueMatchSequencesnoncasc = cell(1,length(coefs));
corSpVec = zeros(1,length(sf));
for i=1:length(coefs) % for each dataset (later move to a function)
    is_distinct_noncasc{i} = zeros(length(sf),size(coefs{i}.allCoefs,2));
    for j=1:length(sf)
        import Core.Discriminative.disc_true;
        [is_distinct_noncasc{i}(j,:), numMatchingSpecies, uniqueMatchSequencesnoncasc{i}{j}] = disc_true(theoriesNoncasc{i}{j}, idSpecies);
        corSpVec(j) = corSpVec(j) + sum(cell2mat(uniqueMatchSequencesnoncasc{i}{j}(find(is_distinct_noncasc{i}(j,:))))~=corSp);
    end
end

%%
ix = 1;
figure; tiledlayout(4,2)
nexttile([1 2])
imagesc(is_distinct_noncasc{ix})
colormap gray
xlabel('Barcode nr.')
% ylabel('Re-scaling factor nr.')
title(' (A) Disc. barcodes for individual re-scaling factor')
nexttile([1 2])

imagesc(is_distinct_casc{ix})
colormap gray
xlabel('Barcode nr.')
% ylabel('Re-scaling factor nr.')
title(' (B) Disc. barcodes for cascading re-scaling factor')
nexttile([2 1])

numDiscBars = zeros(1,length(is_distinct_casc));
numDiscBarsAll = zeros(1,length(is_distinct_casc));
numDiscBarsInd = zeros(1,length(is_distinct_casc));

for ix=1:length(is_distinct_casc)
    sz = size(is_distinct_casc{1},2);
    numDiscBars(ix) = sum(is_distinct_casc{ix}(end,:))/sz;
    numDiscBarsAll(ix) =  sum(sum(is_distinct_casc{ix})>0)/sz;
    numDiscBarsInd(ix) =  sum(sum(is_distinct_noncasc{ix})>0)/sz;
end

% figure
boxplot([numDiscBars' numDiscBarsAll' numDiscBarsInd'],'Labels',{'Disc.','Disc. casc','Disc ind.'})
title(' (C) Percentage of discriminatives, different approaches')


%%


%% cascading
%%

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



%%
[is_distinct_noncascST,is_distinct_cascST,corSpVecST,corSpVecCascST] = Core.casc_stouffer(scoresStouffer,sfCascade,idSpecies,0.6,sf);

ix = 1;
figure; tiledlayout(4,2)
nexttile([1 2])
imagesc(is_distinct_noncascST{ix})
colormap gray
xlabel('Barcode nr.')
% ylabel('Re-scaling factor nr.')
title(' (A) Disc. barcodes for individual re-scaling factor')
nexttile([1 2])

imagesc(is_distinct_cascST{ix})
colormap gray
xlabel('Barcode nr.')
% ylabel('Re-scaling factor nr.')
title(' (B) Disc. barcodes for cascading re-scaling factor')
nexttile([2 1])

numDiscBars = zeros(1,length(is_distinct_cascST));
numDiscBarsAll = zeros(1,length(is_distinct_cascST));
numDiscBarsInd = zeros(1,length(is_distinct_cascST));

for ix=1:length(is_distinct_cascST)
    sz = size(is_distinct_cascST{1},2);
    numDiscBars(ix) = sum(is_distinct_cascST{ix}(end,:))/sz;
    numDiscBarsAll(ix) =  sum(sum(is_distinct_cascST{ix})>0)/sz;
    numDiscBarsInd(ix) =  sum(sum(is_distinct_noncascST{ix})>0)/sz;
end

% figure
boxplot([numDiscBars' numDiscBarsAll' numDiscBarsInd'],'Labels',{'Disc.','Disc. casc','Disc ind.'})
title(' (C) Percentage of discriminatives, different approaches')




% elts = zeros(1,length(rezMaxMP)*length(rezMaxMP{1})*length(rezMaxMP{1}{1}));
% i=1;
% for idx=1:length(rezMaxMP)
%     for idy=1:length(rezMaxMP{idx})
%         for idz = 1:length(rezMaxMP{idx}{idy})
%             elts(i) = rezMaxMP{idx}{idy}{idz}.maxcoef(1)< 0.5;
%             i=i+1;
%         end
%     end



% cellfun(@(z) cellfun(@(y) cellfun(@(x) x{1}.maxcoef(1),y),z),rezMaxMP)
tic
    % allCoefsSingle = cellfun(@(barix) cellfun(@(y) cellfun(@(x) x.maxcoef(1),y),barix,'un',false),rezMaxMP,'un',false);
toc

% allCoefsSingle = arrayfun(@(x) arrayfun(@(y) arrayfun(@(z) rezMaxMP{z}{y}{x}.maxcoef(1) ,1:length(rezMaxMP)),1:length(rezMaxMP{1}),'un',false),1:length(rezMaxMP{1}{1}),'un',false);




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
