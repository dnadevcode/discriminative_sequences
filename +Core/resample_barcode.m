function [outputArg1,outputArg2] = resample_barcode(rezMax,refNums, barGen, bestBarStretch, theoryStruct)
    % resample barcode to get variance & std..
    
    goodbars = bG{idxRun};
    goodbars = goodbars(51);
    w = 300;
    [rezMax,bestBarStretch,bestLength,discSpecies,theoryStruct] = local_alignment_assembly(goodbars, nmbp,w);

    maxCoefs = arrayfun(@(x) rezMax{x}{1}.maxcoef(1),cell2mat(refNums(:)));

    [sortMax,sortMaxId] = sort(maxCoefs,'desc');

    idx = 1;
    idx1 = 1;
    quick_visual_plot(1,refNums{1}(idx),barGen,rezMax,bestBarStretch,theoryStruct)
    
    thrIdx = refNums{1}(idx);
    
    lenBarTested = length(barGen{idx1}.rawBarcode);
    bar1 = interp1(barGen{idx1}.rawBarcode, linspace(1,lenBarTested,lenBarTested*bestBarStretch{thrIdx}(idx1)));
    barB = barGen{idx1}.rawBitmask(round(linspace(1,lenBarTested,lenBarTested*bestBarStretch{thrIdx}(idx1))));
    
    
    % bar1 = imresize(barGen{idx1}.rawBarcode(barGen{idx1}.rawBitmask),'Scale' ,[1 bestBarStretch{thrIdx}(idx1)]) ;
    if rezMax{thrIdx}{idx1}.or(1) == 2
        bar1 = fliplr(bar1);
    end
    bar1 = bar1(barB);
    
    {theoryStruct([refNums{1}(sortMaxId)]).name}'
    
    thr = theoryStruct(thrIdx).theoryBarcode;
    
    pos = find(barGen{idx1}.rawBitmask==1,1,'first');
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
    % s = RandStream('mlfg6331_64'); 
    
    NN = 1000;
    pccScore = zeros(1,NN);
    for i=1:NN
        y = datasample(s,1:numWindows,numWindows,'Replace',true);
        b1 = bar1(R(:,y));
        b2 = bar2(R(:,y));
        pccScore(i) = zscore(b1(:)')*zscore(b2(:))/length(b1(:));
    
    end
    % 
    % figure,plot(zscore(b1(:)));
    % hold on
    % plot(zscore(b2(:)))
    
    mean(pccScore)
    3*std(pccScore)
    
    figure,histogram(pccScore)
    xlabel('CC')