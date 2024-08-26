[hcatheory.sets,hcatheory.names] = Core.Default.read_default_sets('hcasets.txt');

hcaSets = hcatheory.sets.default;

hcaSets.folder = {'C:\Users\Lenovo\postdoc\DATA\Mapping\FASTAS\DA32087.fasta' } ;

hcaSets.computeBitmask = 0;
hcaSets.meanBpExtNm = 0.34;



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


 yoyoConst = 26;
 netropsinConst = 0.4;

yoyo = 25;
netrospin = 1;
hcaSets.theoryGen.method = 'custom';
hcaSets.computeFreeConcentrations = 1;

% hcaSets.pattern = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\new_binding_constants\netroconstants5.txt';
hcaSets.pattern

 % A => 1 , C => 2, G => 3,  T(U) => 4

import CBT.Hca.Core.Theory.choose_cb_model;
[hcaSets.model ] = choose_cb_model(hcaSets.theoryGen.method,hcaSets.pattern, yoyoConst, netropsinConst);

%%


N = 4;

allSeq = arrayfun(@(x)  dec2base(x,4,N), 0:4^N-1, 'UniformOutput', false);
allDna = cellfun(@(x) arrayfun(@(y) str2double(x(y)),1:N)+1, allSeq,'UniformOutput',false);

countATs = cellfun(@(x) sum(x==1)+sum(x==4),allDna);
[sortv,sortord] = sort(countATs,'descend');
sortedSubseq = allDna(sortord);
sortedDNASubseq = cellfun(@(x) int2nt(x),sortedSubseq,'UniformOutput',false)';

sortedVals = cellfun(@(x) hcaSets.model.netropsinBindingConstant(x(1),x(2),x(3),x(4)),sortedSubseq);

% sortedVals = cellfun(@(x) hcaSets.model.netropsinBindingConstant(x(1),x(2),x(3),x(4),x(5)),sortedSubseq);
% vec = arrayfun(@(x1) )
% 
% v1 = 
% % out = tensorprod(eye(4,4),eye(4,4)')

[sortOp,sortid] =  sort(hcaSets.model.netropsinBindingConstant(:),'descend');
[ids ] = arrayfun(@(x) ind2sub([4 4 4 4],x),sortid);

% is this sorted correctly? number of ats
figure,plot(countATs(ids))
figure,plot(sortedVals)


%%
nucleotides = 1:4;
allSeqs = arrayfun(@(x) x, nucleotides )

% theories names
theories = hcaSets.folder;

hcaSets.resultsDir = fullfile(fileparts(theories{1}),'theoryOutput'); % we always save theories output in the same folder as provided data


% make theoryData folder
[~,~] = mkdir(hcaSets.resultsDir);

% compute free concentrations
import CBT.Hca.Core.Theory.compute_free_conc;
hcaSets = compute_free_conc(hcaSets);


%%
% Now sort ATCG's
% it = 1;
% N = 4;
% vec = zeros(256,4);
% for i=1:N
%     for j=1:N
%         for l=1:N
%             for n=1:N
%                 vec(it,:)=[i j l n];
%                 it = it + 1;
%             end
%         end
%     end
% end
