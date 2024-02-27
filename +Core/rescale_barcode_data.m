function [barGen] = rescale_barcode_data(barGen,stretchFactors)

%   Args:
%       barGen - barcode cell structure
%       stretchFactors - stretch factors
%   Returns:
%       barGen - with included rescaled struct

for i=1:length(barGen)
    lenBarTested = length(barGen{i}.rawBarcode);

    for j=1:length(stretchFactors)
        barGen{i}.rescaled{j}.rawBarcode =  interp1(barGen{i}.rawBarcode, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
        barGen{i}.rescaled{j}.rawBitmask  = barGen{i}.rawBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));
    end

end

end

