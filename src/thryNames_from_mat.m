function [uniqueSpeciesNames,idSpecies,thryNames,theoryStruct] = thryNames_from_mat(matFile)
if nargin < 1
    matFile = '/export/scratch/albertas/download_dump/single/theoryOutput/theoryGen_0.34_110_300_0_2024-04-24_19_10_40_session.mat';
end
load(matFile,'theoryNames');

thryNames = theoryNames;

import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name(thryNames);

if nargout >=4
    sets.theoryFile{1} = matFile;
    sets.theoryFileFold{1} = '';
    sets.theory.precision = 5;
    sets.theory.theoryDontSaveTxts = 1;
    import CBT.Hca.UI.Helper.load_theory;
    theoryStruct = load_theory(sets);
end

end

