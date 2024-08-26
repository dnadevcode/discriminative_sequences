function [hcaSets,timestamp] = set_def(folder,nmbp)

[hcatheory.sets,hcatheory.names] = Core.Default.read_default_sets('hcasets.txt');

hcaSets = hcatheory.sets.default;

if nargin < 1
    hcaSets.folder = {'C:\Users\Lenovo\postdoc\DATA\Mapping\FASTAS\DA32087.fasta' } ;
else
    hcaSets.folder = folder;
end

hcaSets.computeBitmask = 0;

if nargin < 2
    hcaSets.meanBpExtNm = 0.34;
else
    hcaSets.meanBpExtNm = nmbp;
end


% make sets compatible with prev structure
hcaSets.theoryGen.method = hcaSets.method;
hcaSets.theoryGen.computeFreeConcentrations =  hcaSets.computeFreeConcentrations;
hcaSets.theoryGen.concN = hcaSets.concN;
hcaSets.theoryGen.concY = hcaSets.concY;
hcaSets.theoryGen.concDNA = hcaSets.concDNA;
hcaSets.theoryGen.psfSigmaWidth_nm = hcaSets.psfSigmaWidthNm;
hcaSets.theoryGen.pixelWidth_nm = hcaSets.pixelWidthNm;
hcaSets.theoryGen.meanBpExt_nm = hcaSets.meanBpExtNm;
hcaSets.theoryGen.isLinearTF = hcaSets.isLinearTF;

hcaSets.lambda.fold  =  hcaSets.fold ;
hcaSets.lambda.name  =  hcaSets.name ;
hcaSets.theoryGen.k =  max(2.^15,2.^(hcaSets.k));
hcaSets.theoryGen.m  = min(2.^15,2.^(hcaSets.m));
hcaSets.theoryGen.computeBitmask = hcaSets.computeBitmask;
% timestamp to add to theories name
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');



end

