function [allSeqA] = find_specific_seq(fastaFold, seqNameToTest, matFile)

    if nargin < 1
        fastaFold = '/export/scratch/albertas/download_dump/single/*.fasta';
    %     seqNameToTest = 'Escherichia coli';
        seqNameToTest = 'Streptococcus pyogenes';
        matFile = '/export/scratch/albertas/download_dump/single/theoryOutput/theoryGen_0.34_110_300_0_2024-04-24_19_10_40_session.mat';

    end

    [uniqueSpeciesNames,idSpecies,thryNames] = thryNames_from_mat(matFile);


    idSeq = find(cellfun(@(x) isequal(x,seqNameToTest),uniqueSpeciesNames)); % need to have species names

    allSeqA = find(idSpecies==idSeq);

    %% if we want to copy these sequences somewhere
    %         sequences = dir(fastaFold);
    
    %     numel(allSeqA)
    %     arrayfun(@(x) copyfile(fullfile(sequences(allSeqA(x)).folder,sequences(allSeqA(x)).name),...
    %         'copyhere'),1:10);
end

