function [outputArg1,outputArg2] = block_bootstrapping(inputArg1,inputArg2)


    %%
thrIdx = refNums{selRef}(idx);%length(theoryStruct)-1; %refNums{1}(idx);
    
    lenBarTested = length(barGenRun{idx1}.rawBarcode);
    bar1 = interp1(barGenRun{idx1}.rawBarcode, linspace(1,lenBarTested,lenBarTested*bestBarStretch{thrIdx}(idx1)));
    barB = barGenRun{idx1}.rawBitmask(round(linspace(1,lenBarTested,lenBarTested*bestBarStretch{thrIdx}(idx1))));
    
    
    % bar1 = imresize(barGen{idx1}.rawBarcode(barGen{idx1}.rawBitmask),'Scale' ,[1 bestBarStretch{thrIdx}(idx1)]) ;
    if rezMax{thrIdx}{idx1}.or(1) == 2
        bar1 = fliplr(bar1);
    end
    bar1 = bar1(barB);
    
    {theoryStruct([refNums{selRef}]).name}'
    
    thr = theoryStruct(thrIdx).theoryBarcode;
    
    pos = find(barGenRun{idx1}.rawBitmask==1,1,'first');
    bar2 = thr(rezMax{thrIdx}{idx1}.pos+pos-1:rezMax{thrIdx}{idx1}.pos+length(bar1)-1+pos-1);
    % 
    zscore(bar1(:)',1)*zscore(bar2(:),1)/length(bar1(:))
    
    figure,plot(zscore(bar1,1));
    hold on
    plot(zscore(bar2))
    
    % now split into ~50px windows
    wWid = 20;
    numWindows = floor(length(bar1)/wWid);
    
    R = reshape(1:floor(length(bar1)/numWindows)*numWindows,floor(length(bar1)/numWindows),[]);
%     s = RandStream('mlfg6331_64'); 
    
    NN = 1000;
    pccScore = zeros(1,NN);
    for i=1:NN
        y = datasample(s,1:numWindows,numWindows,'Replace',true);
        b1 = bar1(R(:,y));
        b2 = bar2(R(:,y));
        pccScore(i) = zscore(b1(:)')*zscore(b2(:))/length(b1(:));
    
    end

end

