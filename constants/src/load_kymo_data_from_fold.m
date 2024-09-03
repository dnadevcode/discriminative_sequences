function [kymoStructs,barN,twoList,bG,expPar,fastaFileF] = load_kymo_data_from_fold(dirName, refsNames,allTheoryFold)
    
    if nargin < 3
        allTheoryFold = '/export/scratch/albertas/download_dump/single/*.fasta';
    end

    import Helper.get_all_folders;
    [barN, twoList] = get_all_folders(dirName);
    
    %% get experiment closest theories
    refsNames = importdata(refsNames);

    import Helper.load_kymo_data;

    
    rez = cell(1,size(twoList,1));
    bG = cell(1,size(twoList,1));
    parfor expNr = 1:size(twoList,1);
        twoList(expNr,:)
        numBars =  barN{twoList(expNr,1)}(twoList(expNr,2))
    
        thry = refsNames{twoList(expNr,1)};
    
        %% get theory location
        dr = dir(allTheoryFold);
    
        thryLoc = find(arrayfun(@(x) ~isempty(strfind(dr(x).name, thry)),1:length(dr)));
    
        fastaFileF{expNr} = fullfile(dr(thryLoc).folder,dr(thryLoc).name);
    
        %% Now load the experiments
        [kymoStructs{expNr}, bG{expNr}, expPar{expNr}.nmpx, expPar{expNr}.nmbp] = load_kymo_data(dirName,1,twoList(expNr,1),twoList(expNr,2));
    
    end

end

