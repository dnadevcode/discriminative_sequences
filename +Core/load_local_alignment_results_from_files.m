function [rM,bnames,mpval,thryNames] = load_local_alignment_results_from_files(foldCalc)
    %   load_local_alignment_results_from_files
    %   
    %   Args:
    %       foldCalc - folder with kymographs
    %
    %   Returns:
    %       rM - result struct
    %       bnames - barcodes in the curren results struct
    %       mpval - lengths of local alignments, 0 - full alignment.

    files = dir(fullfile(foldCalc,'..','*.txt'));
    
    
    N = length(files);
    
    mpval = zeros(1,N);
    for i=1:N
        f1= strsplit(files(i).name,'w=');
        if length(f1) ~= 1
            f2 = strsplit(f1{2},'_');
            mpval(i) = str2double(f2{1});
        end
    end
    
    rM = cell(1,N);
    bnames = cell(1,N);
    
    % in case output structure changes (i.e. different number of outputs,
    % format, change load_coefs accordingly).
    for i=1:N
        [rM{i},bnames{i},thryNames{i}] = Core.load_coefs(fullfile(files(i).folder,files(i).name));
    end

end

