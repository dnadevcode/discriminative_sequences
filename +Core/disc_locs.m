function [refNums,allNums] = disc_locs(rezMax)

refNums = cell(1,length(rezMax{1}));
for barIx = 1:length(rezMax{1});
%     toc
    allCCs = cellfun(@(x) x{barIx}.maxcoef(1),rezMax);

% figure,plot(allCCs)
% xlabel('Theory bar')
% ylabel('Score')

    [a,sortedid] = sort(allCCs,'desc','MissingPlacement','last');
    maxs = a(1);
    cdiff = 0.05;
    selectedRef = sortedid(find(a>a(1)-cdiff));
    refNums{barIx} = selectedRef;
end

allNums = cellfun(@(x) length(x),refNums);
end

