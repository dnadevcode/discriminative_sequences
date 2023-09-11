function [pccScore] = block_bootstrapping(b1,bT,overlapStruct,k,iy, wLocal, wBoostrapping,s)
    % block_bootstrapping including mp
    %
    %   Args:
    %       b1 - experimental barcode
    %       bT - theory barcode
    %       overlapStruct - the results structure
    %       k - idx first barcode, usually 1
    %       iy - idx second barcode, usually 2
    %       wLocal - local window width. 0 - global
    %       wBootstrapping - boostrapping window, here set to 20
    %       s - random stream
    %   Returns:
    %       pccScore : vector of bootstrapped values

    %
    if nargin < 7
        wBoostrapping = 20;   % block window width
        s = RandStream('mlfg6331_64'); 
    end

    % convert to matching structures
     b1Struct = cell2struct([cellfun(@(x) double(x.rawBarcode),b1,'un',false);...
        cellfun(@(x) x.rawBitmask,b1,'un',false)]',{'rawBarcode','rawBitmask'},2);

    barcodeGenTStruct = struct('rawBarcode',bT(1).rawBarcode,'rawBitmask',true(1,length(bT(1).rawBarcode)));

     barStruct = [b1Struct barcodeGenTStruct];

     if wLocal == 0 % todo: move to get_local_subbarcode_pair
        lenBarTested = length(barStruct(1).rawBarcode);
        bar1 = interp1(barStruct(1).rawBarcode, linspace(1,lenBarTested,lenBarTested*overlapStruct.bestBarStretch{1}));
        barB = barStruct(1).rawBitmask(round(linspace(1,lenBarTested,lenBarTested*overlapStruct.bestBarStretch{1})));
        
        if overlapStruct.rezMax{1}{1}.or == 2
            bar1 = fliplr(bar1);
        end
        bar1 = bar1(barB);
            
        pos = overlapStruct.rezMax{1}{1}.secondPos;
        bar2 = barStruct(2).rawBarcode(overlapStruct.rezMax{1}{1}.pos(1)+pos-1:overlapStruct.rezMax{1}{1}.pos(1)+length(bar1)-1+pos-1);
    % 
     else
    
        [bar1, bar2] = get_local_subbarcode_pair(barStruct,overlapStruct,k,iy);

     end

%     zscore(bar1(:)',1)*zscore(bar2(:),1)/length(bar2(:))
%     
% figure,plot(zscore(bar1,1));
% hold on
% plot(zscore(bar2,1))
%     
    
    % now split into ~50px windows. skips the last window if not divisible
    numWindows = floor(length(bar1)/wBoostrapping);
    
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

