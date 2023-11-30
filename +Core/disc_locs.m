function [refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locs(rezMax,stouffer)
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

    if nargin < 2
        stouffer = 0;
    end
    cdiff = 0.05;

    numTheories = length(rezMax{1});

    
    refNums = cell(1,numTheories);
    bestCoefs =  cell(1,numTheories);
    refNumBad =  zeros(1,numTheories);
    bestCoefsBad =  zeros(1,numTheories);

    for barIx = 1:length(rezMax{1})
        allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMax);

        if stouffer
            allCCStouffer = cellfun(@(x) x{barIx}.stoufferScore(1),rezMax); 
            [a,sortedid] = sort(allCCStouffer,'desc','MissingPlacement','last');
            allCCsS = allCCs(sortedid);
            %         maxs = a(1);
            
            discLocations = (find(allCCsS>allCCsS(1)-cdiff)); % same cdiff until we update
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

        else
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
    
    % figure,plot(allCCs)
    % xlabel('Theory bar')
    % ylabel('Score')
    


    end
    
    allNums = cellfun(@(x) length(x),refNums);
end

