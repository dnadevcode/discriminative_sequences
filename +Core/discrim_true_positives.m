function [truePositives,discSpecies,discAll,allNums,refNums,signMatch,fp,positives] = discrim_true_positives(rezMax,speciesLevel, idc)

%Returns:   
%   truePositives - true positives based on allspecies
%   discSpecies - how many strains is it discriminative to
%   discAll - binary vector of all strains its discriminative to
%   allNums - all discriminative locations
%   refNums - refnums of discriminative locations
%   signMatch - numver of unique discriminative locations
import Core.disc_locs;
[refNums,allNums] = disc_locs(rezMax);

signMatch = find(logical(allNums ==1));
% refNums(signMatch)
% theoryStruct([cell2mat(refNums(signMatch))]).name;


allSpecies = find(speciesLevel);
allSpecUnique = unique(idc(find(speciesLevel)));
try
idc(idc==allSpecUnique(2)) = allSpecUnique(1);
catch
end
discAll = cellfun(@(x) ismember(x,allSpecies),refNums,'UniformOutput',false);

discSpecies = cellfun(@(x) sum(ismember(x,allSpecies)==0),refNums,'UniformOutput',true);

truePositives = sum(discSpecies==0); % also count false positives?;

discAny = cellfun(@(x) unique(idc(x)),refNums,'UniformOutput',false);

positives = cellfun(@(x) length(x),discAny);
fp = sum(positives==1)-truePositives;

end

