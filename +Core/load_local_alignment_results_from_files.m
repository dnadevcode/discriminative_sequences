function [rM,bnames,mpval,thryNames,files] = load_local_alignment_results_from_files(foldCalc,file)
    %   load_local_alignment_results_from_files
    %   
    %   Args:
    %       foldCalc - folder with kymographs
    %
    %   Returns:
    %       rM - result struct
    %       bnames - barcodes in the curren results struct
    %       mpval - lengths of local alignments, 0 - full alignment.

    if nargin < 2 || file == 0
        files = dir(fullfile(foldCalc,'..','*.txt'));
    else
        files = dir(file);
    end
    
    
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
    thryNames = cell(1,N);

%     ds = datastore(fullfile(files(i).folder,files(i).name));
%     tt = tall(ds);
%     thryNames{1} = gather(tt.Var1);   
%     bnames{1} = ds.SelectedVariableNames(2:4:end);
%     
    
    % in case output structure changes (i.e. different number of outputs,
    % format, change load_coefs accordingly).
    for i=1:N % could be parfor, but probably no speedup
        [rM{i},bnames{i},thryNames{i}] = Core.load_coefs(fullfile(files(i).folder,files(i).name));
    end

end

