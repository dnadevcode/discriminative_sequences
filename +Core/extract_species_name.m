function [speciesLevel] = extract_species_name(theoryStruct)
% names = [theoryStruct(1:10).name];
names = arrayfun(@(x) theoryStruct(x).name,1:length(theoryStruct),'un',false);

ecoliTheories = cellfun(@(x) ~isempty(strfind(x,'Escherichia coli')),names);
shigellaTheories = cellfun(@(x) ~isempty(strfind(x,'Shigella')),names);

speciesLevel = ecoliTheories+shigellaTheories;
% 
% klepTheories = cellfun(@(x) ~isempty(strfind(x,'Klebsiella pneumoniae')),names);
% 
% 
% 
% numEcolis = sum(ecoliTheories);
% numShigella = sum(shigellaTheories);


end

