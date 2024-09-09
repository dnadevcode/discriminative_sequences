function [distinctoutput, mostCommonSequence,mostCommonRep,namesDir,thryDiscName] = ...
    get_best_theory(savedir, refNums, idSpecies,barN, twoList,foldname,thryNames, sepname, sfLevel, cdiff)
    % get id for best theory for each of the experiment folders

    if nargin < 10
        sfLevel = 15; 
        cdiff = 0.05;
    end

    if nargin < 8
        sepname = 'oldhca';
        split2 = 'Sample';
    end


    dirs = dir(fullfile(savedir,'*.mat'));
    % distinctoutput
   lf = @(x,y) load(fullfile(x(y).folder,x(y).name));

   mostCommonSequence = zeros(1,length(dirs));
   mostCommonRep = zeros(1,length(dirs));
   namesDir= cell(1,length(dirs));
    thryDiscName= cell(1,length(dirs));
    distinctoutput = cell(1,length(dirs));

%%
    for nr = 1:length(dirs);
        nr

        maxData = lf(dirs,nr);
        num_disc = zeros(floor(size(maxData.matAllCoefs,2)/2+1),length(cdiff));
        numFP = zeros(floor(size(maxData.matAllCoefs,2)/2+1),length(cdiff));
        
        maxCoef = cell(1,size(maxData.matAllCoefs,1));
        for barid =1:size(maxData.matAllCoefs,1)
            [singleCoef , singlePos ] =  max(maxData.matAllCoefs(barid,sfLevel:end-sfLevel+1,:),[],2);
            pos  = squeeze(singlePos)';
            maxCoef{barid} =  squeeze(singleCoef);
          
        end
        
        import Core.Discriminative.disc_locations;
        [refNums, allNums, bestCoefs,refNumBad,bestCoefsBad] = disc_locations({maxCoef}, cdiff);
        
        import Core.Discriminative.disc_true;
        [is_distinct,numMatchingSpecies,uniqueMatchSequences] = disc_true(refNums, idSpecies);
        
        % mostCommonSequence(nr) = mode(arrayfun(@(x) refNums{x}(1),find(is_distinct));
        % mostCommonRep(nr) = sum(arrayfun(@(x) refNums{x}(1),find(is_distinct))==mostCommonSequence(nr));

        
        mostCommonSequence(nr) = mode(cell2mat(arrayfun(@(x) refNums{x},find(is_distinct),'un',false)')); % todo : weight by coefficien/position?
        mostCommonRep(nr) = sum(cell2mat(arrayfun(@(x) refNums{x},find(is_distinct),'un',false)')==mostCommonSequence(nr));


        % %%
        % dirn = strsplit(dirs(nr).name,'oldhca');
        % dirn = str2num(dirn{1});
        % nameFold = foldname{twoList(dirn,1)}{twoList(dirn,2)};
        % thryNames{mostCommonSequence(nr)}

        if mostCommonRep(nr)>0

            dirn = strsplit(dirs(nr).name,['' ...
                  '' ...
                  sepname]);
            dirn = str2num(dirn{1});
            nameFold = foldname{twoList(dirn,1)}{twoList(dirn,2)};
            nameFold = strsplit(nameFold,'_');
            nameFold = nameFold{2};
            % nameFold = strsplit(nameFold,split2);
            % nameFold = strsplit(nameFold{2},'-');
            % nameFold = str2num(nameFold{1})
    
            namesDir{nr} = nameFold;
    
            thryDisc = strsplit(thryNames{mostCommonSequence(nr)},' ');
            thryDiscName{nr} = thryDisc{1}(2:end);

        end


    end


    %% 
    [namesDir' thryDiscName' arrayfun(@(x) x,mostCommonRep,'un',false)']


end

