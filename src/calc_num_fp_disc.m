function [num_disc, numFP, numUnique] = calc_num_fp_disc(lf,dirs,nr,cdiff,idSeq,idSpecies)
% calculate number of discriminative barcodes and false positives

maxData = lf(dirs,nr);


% sF = maxData.sets.theory.stretchFactors;
% 
% cdiff = 0.01:0.01:0.1;
% idSeq = 4717;



num_disc = zeros(floor(size(maxData.matAllCoefs,2)/2+1),length(cdiff));
numFP = zeros(floor(size(maxData.matAllCoefs,2)/2+1),length(cdiff));
numUnique= zeros(floor(size(maxData.matAllCoefs,2)/2+1),length(cdiff));

for j=1:floor(size(maxData.matAllCoefs,2)/2+1)
    sfLevel = j;
    maxCoef = cell(1,size(maxData.matAllCoefs,1));
    for barid =1:size(maxData.matAllCoefs,1)
        [singleCoef , singlePos ] =  max(maxData.matAllCoefs(barid,sfLevel:end-sfLevel+1,:),[],2);
        pos  = squeeze(singlePos)';
        maxCoef{barid} =  squeeze(singleCoef);
      
    end 
    % load([fold,'maxcoefs',num2str(j),'.mat']); % if presaved
    
    for i=1:length(cdiff)
        import Core.Discriminative.disc_locations;
        [refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoef}, cdiff(i));
        
        import Core.Discriminative.disc_true;
        [is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);
        
        num_disc(j,i) = sum(is_distinct);
        
        numUnique(j,i) =sum(cellfun(@(x) length(x),refNums)==1);
        
        numFP(j,i) = sum(idSpecies(arrayfun(@(x) refNums{x}(1),find(is_distinct))) ~= idSeq);
    
    end

end

end

