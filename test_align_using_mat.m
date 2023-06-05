
thryFile = '/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files May 2022/Ref_Theory_220502_130nmPixel_0.34nmPERbp_HCAv4.4.mat';

% load(thryFile);

tic
sets.w = 300;
sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
sets.theory.theoryDontSaveTxts = 1;
sets.theoryFile{1} = thryFile;
% compare theory to experiment

import CBT.Hca.UI.Helper.load_theory;
theoryStruct = load_theory(sets);

sets.theory.nmbp = 0.243;

tic
import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct, sets );
toc
% create struct
% theoryStruct = cell2struct([theoryGen.theoryBarcodes;...
%     theoryGen.theoryBitmasks;arrayfun(@(x) 0,1:length(theoryGen.theoryBarcodes),'un',false)]',{'theoryBarcodes','theoryBitmasks', 'isLinear'},2);

% barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barGen,'un',false);...
%     cellfun(@(x) x.rawBitmaskbarGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

import CBT.Hca.Core.Comparison.compare_distance;
[rezMax,bestBarStretch,bestLength] = compare_distance(barGen,theoryStruct, sets, [] );