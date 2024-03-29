function [] = export_coefs_stouffer_norm(theoryStruct,rezMax, barGen,matDirpath)

    % export_cc_vals_table
    % only for single barcodes, dont include the consensus here
    
%     fasta = cell(1,length(theoryStruct));
%     for i =1:length(theoryStruct)
%         fasta{i} = theoryStruct(i).name;
%     end
    %fasta = cellfun(@(x) strrep(x.filename,'.txt',''), theoryStruct,'UniformOutput',false);
    thrLen = arrayfun(@(x) theoryStruct(x).length, 1:length(theoryStruct));
    fasta = arrayfun(@(x) theoryStruct(x).name, 1:length(theoryStruct),'un',false);

    T = table(fasta');

    for i=1:length(barGen)
        maxccoef = cell2mat(cellfun(@(x) x{i}.stoufferScoreNorm(1), rezMax,'UniformOutput',0));
        lengthPx = cell2mat(cellfun(@(x) x{i}.lengthMatch,rezMax,'UniformOutput',0));
        pos = cell2mat(cellfun(@(x) x{i}.pos(1), rezMax,'UniformOutput',0));
        stretch =  cell2mat(cellfun(@(x) x{i}.bestBarStretch,rezMax,'UniformOutput',0));
        for j=1:length(pos)
           if pos(j)<= 0
               pos(j) = pos(j)+thrLen(j);
           end
        end
        maxcc = maxccoef;
        
        [d,name,ext] = fileparts(barGen{i}.name);

        N = matlab.lang.makeValidName(name);
        
        T2 = table(maxcc',lengthPx', pos',stretch' ,'VariableNames',{N,strcat(['len_'  num2str(i)]),strcat(['pos_'  num2str(i)]),strcat(['stretch_'  num2str(i)])});
        T = [T T2];
    end
    disp('Saving ccvals table');
    import CBT.Hca.Export.export_cc;
    export_cc(T, matDirpath);      


end

