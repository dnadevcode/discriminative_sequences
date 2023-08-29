function [pccScore] = block_bootstrapping(barStruct,overlapStruct,k,iy,wWid,s)
    % block_bootstrapping MP
    if nargin < 5
        wWid = 20;   % block window width
        s = RandStream('mlfg6331_64'); 
    end


    [bar1, bar2] = get_local_subbarcode_pair(barStruct,overlapStruct,k,iy);
%     
% figure,plot(zscore(subBar1,1));
% hold on
% plot(zscore(subBar2))
%     
    
    % now split into ~50px windows
    numWindows = floor(length(bar1)/wWid);
    
    R = reshape(1:floor(length(bar1)/numWindows)*numWindows,floor(length(bar1)/numWindows),[]);
     
    
    NN = 1000;
    pccScore = zeros(1,NN);
    for i=1:NN
        y = datasample(s,1:numWindows,numWindows,'Replace',true);
        b1 = bar1(R(:,y));
        b2 = bar2(R(:,y));
        pccScore(i) = zscore(b1(:)')*zscore(b2(:))/length(b1(:));
    
    end

end

