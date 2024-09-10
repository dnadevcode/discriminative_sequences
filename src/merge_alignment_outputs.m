function [] = merge_alignment_outputs(saveDir,newDir, twoList)

%     dirs = dir(fullfile(saveDir,'*.mat'));
    % distinctoutput
   lf = @(x,y) load(fullfile(x(y).folder,x(y).name),'matAllCoefs');
% 
%    mostCommonSequence = zeros(1,length(dirs));
%    mostCommonRep = zeros(1,length(dirs));
%    namesDir= cell(1,length(dirs));
%     thryDiscName= cell(1,length(dirs));
%     distinctoutput = cell(1,length(dirs));

    for i=1:max(twoList(:,1))
        curIds = find(twoList(:,1)==i);
    
        dirs = arrayfun(@(x) dir(fullfile(saveDir,[num2str(x),'_','*.mat'])),curIds,'un',true);

        maxData = arrayfun(@(x) lf(dirs,x),1:length(dirs),'un',false);
        data = cellfun(@(x) x.matAllCoefs,maxData,'UniformOutput',false);
        matAllCoefs = cat(1,data{:});
        save([newDir, num2str(i),'_','sf_allcoefs','.mat'],'matAllCoefs','-v7.3');
    end


end

