function [subBar1, subBar2] = get_local_subbarcode_pair(barStruct,overlapStruct,k,iy)
    
    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end

 
    if isfield(overlapStruct,'pA')
        pA = overlapStruct(k,iy).pA ;
        pB = overlapStruct(k,iy).pB ;
        h = overlapStruct(k,iy).h;
        orr = overlapStruct(k,iy).or ;
        bS =  overlapStruct(k,iy).bestBarStretch;
        score = overlapStruct(k,iy).score;
        fullscore = overlapStruct(k,iy).fullscore;
    
       % resizing depends on how barcodes were defined. For PCC, we resize
        % first, and then apply bitmask. For MP?
        bBar = imresize(barStruct(k).rawBarcode(barStruct(k).rawBitmask),'Scale' ,[1 bS]);

    else
        % in this case we have slightly different structure..

        pA = overlapStruct.rezMax{iy-1}{k}.secondPos(1);
        pB = overlapStruct.rezMax{iy-1}{k}.pos(1);
        h = overlapStruct.rezMax{iy-1}{k}.lengthMatch;
        orr = overlapStruct.rezMax{iy-1}{k}.or(1);
        bS = overlapStruct.bestBarStretch{iy-1};
        score = overlapStruct.rezMax{iy-1}{k}.maxcoef(1);
        fullscore = nan; % using matlab's version, we don't calculate this within function

        lenBarTested = length(barStruct(k).rawBarcode);
        bBar = interp1(barStruct(k).rawBarcode, linspace(1,lenBarTested,lenBarTested*bS));
        bBit = barStruct(k).rawBitmask(round(linspace(1,lenBarTested,lenBarTested*bS)));
        bBar(~bBit) = nan;

    end

    aBar = barStruct(iy).rawBarcode(barStruct(iy).rawBitmask);


    if orr~=1
        bBar = fliplr(bBar);
    end


    subBar1 = bBar(pA:pA+h-1);
    subBar2 = aBar(pB:pB+h-1);



end

