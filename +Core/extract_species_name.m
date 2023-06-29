function [speciesLevel,idc] = extract_species_name(theoryStruct)


% names = [theoryStruct(1:10).name];
names = arrayfun(@(x) theoryStruct(x).name,1:length(theoryStruct),'un',false);

ecoliTheories = cellfun(@(x) ~isempty(strfind(x,'Escherichia')),names);
shigellaTheories = cellfun(@(x) ~isempty(strfind(x,'Shigella')),names);

speciesLevel = ecoliTheories+shigellaTheories;

if nargout > 1
    % also return unique identifiers for each species
    species = arrayfun(@(x) strsplit(theoryStruct(x).name,' '),1:length(theoryStruct),'un',false);
    species = cellfun(@(x) x{2},species,'un',false);
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

