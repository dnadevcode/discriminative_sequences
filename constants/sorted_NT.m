function [sortedSubseq,sortv,countATs,orderSeq,sortord,allDna] = sorted_NT(N)


% N = 5;
allSeq = arrayfun(@(x)  dec2base(x,4,N), 0:4^N-1, 'UniformOutput', false);
allDna = cellfun(@(x) arrayfun(@(y) str2double(x(y)),1:N)+1, allSeq,'UniformOutput',false);

countATs = cellfun(@(x) sum(x==1)+sum(x==4),allDna);
[sortv,sortord] = sort(countATs,'descend');
sortedSubseq = allDna(sortord);

orderSeq = zeros(1,numel(allDna));
for i=1:numel(allDna)
    cellInd = num2cell(allDna{i});
    orderSeq(i) = sub2ind(ones(1,N)*4, cellInd{:} );
end


end
