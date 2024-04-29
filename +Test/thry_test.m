% Converting to specific nm/bp ratio

% Check that bytes send and bytes received match with the output theory
% array
tic
ticBytes(gcp);

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct2 = convert_nm_ratio(0.15, theoryStruct, sets );

tocBytes(gcp);
toc

whos('theoryStruct')
1715738244/1024^3
% theoryStruct2.rawBarcode


%  s = builtin('struct', theoryStruct2); 
% whos('s')


% size of the theoryStruct array
onlyBar = {theoryStruct2(:).rawBarcode};

arraysize = sum(cellfun(@(x) length(x),onlyBar))*8;
fprintf("The array will need %.0f bytes in memory", arraysize)

inKB = arraysize/1024;
inMB = arraysize/1024^2;
inGB = arraysize/1024^3;
fprintf("The array will need %.5f KB or %.5f MB or %.5f GB of memory", inKB, inMB, inGB)


%% BG

% rescale. Could remove unnecessary fields here?
import Core.rescale_barcode_data;
[barGen2] = rescale_barcode_data(barGen,sets.theory.stretchFactors);
whos('barGen')
whos('barGen2')

272037/1024^3


%% rezmax
whos('rezMaxMP')
47464704/1024^3*50


numel(allCoefs(:))/1024^2*50