
[hcaSets,timestamp] = set_def();

yoyoConst = 25;
netropsinConst = 0.4;

yoyo = 25;
netrospin = 1;
hcaSets.theoryGen.method = 'custom';
hcaSets.computeFreeConcentrations = 1;
hcaSets.pattern = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\new_binding_constants\netroconstants5.txt';


import CBT.Hca.Core.Theory.choose_cb_model;
[hcaSets.model ] = choose_cb_model(hcaSets.theoryGen.method,hcaSets.pattern, yoyoConst, netropsinConst);

N = 5;

allSeq = arrayfun(@(x)  dec2base(x,4,N), 0:4^N-1, 'UniformOutput', false);
allDna = cellfun(@(x) arrayfun(@(y) str2double(x(y)),1:N)+1, allSeq,'UniformOutput',false);

countATs = cellfun(@(x) sum(x==1)+sum(x==4),allDna);
[sortv,sortord] = sort(countATs,'descend');
sortedSubseq = allDna(sortord);
sortedDNASubseq = cellfun(@(x) int2nt(x),sortedSubseq,'UniformOutput',false)';

sortedVals = cellfun(@(x) hcaSets.model.netropsinBindingConstant(x(1),x(2),x(3),x(4)),sortedSubseq);

%% sort
% sigma = 20;
% constFun = arrayfun(@(x) 30*exp(-x/sigma),1:length(sortedSubseq));
% figure,plot(constFun);
% 
% for i=1:length(sortedSubseq)
%     hcaSets.model.netropsinBindingConstant(sortedSubseq{i}(1),sortedSubseq{i}(2),sortedSubseq{i}(3),sortedSubseq{i}(4)) = constFun(i);
% end

[sortOp,sortid] =  sort(hcaSets.model.netropsinBindingConstant(:),'descend');
[ids ] = arrayfun(@(x) ind2sub([4 4 4 4 4],x),sortid);

% is this sorted correctly? number of ats
figure,plot(countATs(ids))
xlabel('5mers sorted by binding constant')
ylabel('AT content')