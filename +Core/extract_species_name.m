function [speciesLevel,idc] = extract_species_name(theoryStruct,namesSpecies,names)

if nargin < 2
    namesSpecies = {'Escherichia','Shigella'};
end

if nargin < 3
% names = [theoryStruct(1:10).name];
    names = arrayfun(@(x) theoryStruct(x).name,1:length(theoryStruct),'un',false);
end

speciesLevel = zeros(1,length(names));
for i=1:length(namesSpecies);
    speciesTheories = cellfun(@(x) ~isempty(strfind(x,namesSpecies{i})),names);
    speciesLevel = speciesLevel + speciesTheories;
end
% speciesLevel = ecoliTheories+shigellaTheories;

if nargout > 1
    % also return unique identifiers for each species
    species = arrayfun(@(x) strsplit(names{x},' '),1:length(names),'un',false);
%     species = cellfun(@(x) x{2},species,'un',false);
    species = cellfun(@(x) [x{2},' ',x{3}],species,'un',false);

    % cellfun(@(x) unique(x))
    [uc, ~, idc] = unique( species ) ;
end
% 
% klepTheories = cellfun(@(x) ~isempty(strfind(x,'Klebsiella pneumoniae')),names);
% 
% 
% 
% numEcolis = sum(ecoliTheories);
% numShigella = sum(shigellaTheories);


end

