function [rezOut] = local_alignment_assembly(theoryStruct, barGen, w, sets)
    % Tailored to work with kymographs from barcode assembly problem to
    % find which are discriminative, so we skip to bargen

    %   Args:
    %       theoryStruct - theory structure
    %       barGen - barcode to calculate alignments for
    %       w - alignment length
    %
    %   Returns:
    %       rezOut - output results structure (which then could be saved to
    %       a txt file

    if isempty(theoryStruct)
        nmbp = 0.25;
        %% theory loading
        
        import Core.load_theory_structure;
        thryFileIdx = 1; % todo: pass directly the theory file here
        [theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);

    end

    if nargin < 4
        % we following "Strain-level bacterial typing directly from patient
        % samples using optical DNA mapping", i.e. 20 timeframes, 10%
        % stretching at 2.5% intervals, nralign alignment and otsu-based
        % edge detection.
        sets.timeFramesNr = 20;
        sets.theory.stretchFactors = 0.9:0.025:1.1; %as per 
        sets.alignMethod = 1;
        sets.edgeDetectionSettings.method = 'Otsu';
        sets.genConsensus = 0;
        sets.filterSettings.filter = 0;
        sets.dirName = 'output';
        sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    end

    %
    if isempty(w)
        w = 0;
        sets.w = 0;
    end

    rezOut = cell(1,length(w));

    
    % speciesLevel
    
    import CBT.Hca.Core.Comparison.compare_distance;
    
    for wIdx = 1:length(w)
        sets.w = w(wIdx);
        passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask)*sets.theory.stretchFactors(end), barGen) >= sets.w); % which barcodes are long enough?
    
        %
        if sets.w == 0
            sets.comparisonMethod = 'mass_pcc';
        else
            sets.comparisonMethod = 'mpnan';
        end
    
        [rezMaxMP,bestBarStretchMP,bestLengthMP] = compare_distance(barGen(passingThreshBars),theoryStruct, sets, [] );
    
        rezOut{wIdx}.rezMax = rezMaxMP;
        rezOut{wIdx}.bestBarStretch = bestBarStretchMP;
        rezOut{wIdx}.bestLength = bestLengthMP;
    
    %     import Core.export_coefs;
    %     export_coefs(theoryStruct,rezMaxMP,bestBarStretchMP,barGen(passingThreshBars),[sets.dirName, '_MP_w=',num2str(sets.w),'_']);
    %     save([sets.dirName, num2str(sets.w),'_rez.mat'],'rezMaxMP','passingThreshBars','sets');
    end
    
    %% Visual evaluation plots.
    
    % quick_visual_plot(16,9242,barGen,rezMax,bestBarStretch,theoryStruct)
    
    %  super_quick_plot(16,barGen,comparisonStruct,theoryStruct)
    % sigmatches = find(allNums ==1)
    % for i=1:length(sigmatches)
    %     quick_visual_plot(sigmatches(i),9242,barGen,rezMax,bestBarStretch,theoryStruct)
    % end
    
    % cell2mat(refNums(sigmatches))
    % refNums(signMatch)
    % theoryStruct([cell2mat(refNums(signMatch))]).name;


end