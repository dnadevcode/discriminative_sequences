function [refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locs(rezMax)
    % disc_locs - simple code to find locations within cdiff of the max. If
    % they are all from the same species/subspecies, then the barcode is
    % discriminative at that level

    %   Args:
    %       rezMax
    %   Returns:
    %       refNums - best scoring references for each
    %       allNums - number of sequences within the limit
    %       bestCoefs - best coefficients
    %       refNumBad - first theory below cdiff
    %       bestCoefsBad - first coefficient below cdiff

    cdiff = 0.05;

    numTheories = length(rezMax{1});

    
    refNums = cell(1,numTheories);
    bestCoefs =  cell(1,numTheories);
    refNumBad =  zeros(1,numTheories);
    bestCoefsBad =  zeros(1,numTheories);

    for barIx = 1:length(rezMax{1})
        allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMax);
    
    % figure,plot(allCCs)
    % xlabel('Theory bar')
    % ylabel('Score')
    
        [a,sortedid] = sort(allCCs,'desc','MissingPlacement','last');
%         maxs = a(1);

        discLocations = (find(a>a(1)-cdiff));
        selectedRef = sortedid(discLocations);
        refNums{barIx} = selectedRef;
        bestCoefs{barIx} = a(discLocations);
        try
            refNumBad(barIx) = sortedid(discLocations(end)+1);
            bestCoefsBad(barIx) = a(discLocations(end)+1);
        catch
            refNumBad(barIx) = nan;
            bestCoefsBad(barIx) = nan;
        end

    end
    
    allNums = cellfun(@(x) length(x),refNums);
end

