function [scores,pccScore] = local_bootstrap_run(barGenRun,rM,bnames,theoryStruct ,mpval,speciesLevel,idc)
% local_bootstrap_run
%   Args:
%       barGenRun, rM, bnames, theoryStruct , mpval, speciesLevel, idc   
%   Returns:
%       scores - mean( pccScore{m}) - 3*std(pccScore{m});
N = length(rM); % num different lengths

scores = nan(N,4);
pccScore = cell(1,N);

for m = 1:N % which length to analyse
    % Step 1 : find the best scoring theory for each barcode  
    import Core.disc_locs;
    [refNums, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rM{m});

    % calcs ref and allnums again..
    import Core.discrim_true_positives;
        [truePositivesMP,discSpeciesMP,discAllMP,allNumsMP,refNumsMP,signMatchMP, fpMP,positives] = ...
        discrim_true_positives(rM{m}, speciesLevel, idc);

    % intermediate step, check if this barcode is within the output struc
    % (depeends on its length)
    [barsPassThresh,locb] = ismember(cellfun(@(x) matlab.lang.makeValidName(strrep(x.name,'.tif','')),barGenRun,'un',false),bnames{m});

    if barsPassThresh~=0 && sum(bestCoefs{locb})~=0
        
        % Step 2 : recalculate (if needed) best alignment against theory. Useful to check
        % that the alignment was correct. Not needed if rezults is saved as mat
        % file
        
        theoryStructSel = theoryStruct(refNums{locb}(1)); % against the same theory, but we don't know which is correct! so we calcualte it for each length separately
        % barGenRun = barGen(1);
        % w = [200:50:sum(barGenRun{1}.rawBitmask)]; % in practice could run all
        [rezOutRecalc] = local_alignment_assembly(theoryStructSel, barGenRun,mpval(m));
    
        % Step 3 : block bootstrapping
        pccScore{m} = block_bootstrapping(barGenRun, theoryStructSel,rezOutRecalc{1},1, 2,mpval(m));
    
        scores(m,:) = [mean( pccScore{m}) 3*std(pccScore{m}) length(refNums) positives(locb)];
    else
    end
end

end

