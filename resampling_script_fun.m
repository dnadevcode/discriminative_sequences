function [] = resampling_script_fun(dirName,ix)
    
    import Core.load_data_fun;
    [kymoStructs,barGen,sets] = load_data_fun(dirName,ix);
    
    curDirKymos = sets.dirName; % current directory with kymographs
    
    nmbp = sets.nmbp;
    
    import Core.load_theory_structure;
    thryFileIdx = 1; % todo: pass directly the theory file here
    [theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);
    
    
    import Core.extract_species_name;
    [speciesLevel,idc] = extract_species_name(theoryStruct);
    %
    % load alignment result to result struct if it was saves as text files
    import Core.load_local_alignment_results_from_files;
    [rM, bnames, mpval] = load_local_alignment_results_from_files(curDirKymos ); 
    
    
    
    
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    
    scores = cell(1,length(barGen));
    
    for ii=1:length(barGen)
        [scores{ii},pccScore] = local_bootstrap_run( barGen(ii),rM,bnames,theoryStruct ,mpval,speciesLevel,idc);
        scores{ii}
        
        % output:
        import Core.export_coefs_resampling;
        T = export_coefs_resampling(scores{ii}, barGen(ii), mpval, [curDirKymos, '_resampling_table'],timestamp);
    
        save( [curDirKymos, '_resampling_table.mat'],'T');
    end


end

