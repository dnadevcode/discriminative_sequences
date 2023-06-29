function [rezMax,barnames] = load_coefs(coefdata)

    data = importdata(coefdata);

    barnames = data.textdata(1,2:4:end);

    % slow
    tic
    rezMax =  cell(1,size(data.data,1));
    for i=1:size(data.data,1)
        for j=1:size(data.data,2)/4
            rezMax{i}{j}.maxcoef = data.data(i,1+4*(j-1));
            rezMax{i}{j}.lengthMatch = data.data(i,2+4*(j-1));
            rezMax{i}{j}.pos = data.data(i,3+4*(j-1));
            rezMax{i}{j}.bestBarStretch = data.data(i,4+4*(j-1));
        end
    end

    toc

    
%     % export_cc_vals_table
%     % only for single barcodes, dont include the consensus here
%     
% %     fasta = cell(1,length(theoryStruct));
% %     for i =1:length(theoryStruct)
% %         fasta{i} = theoryStruct(i).name;
% %     end
%     %fasta = cellfun(@(x) strrep(x.filename,'.txt',''), theoryStruct,'UniformOutput',false);
%     thrLen = arrayfun(@(x) theoryStruct(x).length, 1:length(theoryStruct));
%     fasta = arrayfun(@(x) theoryStruct(x).name, 1:length(theoryStruct),'un',false);
% 
%     T = table(fasta');
% 
%     for i=1:length(barGen)
%         maxccoef = cell2mat(cellfun(@(x) x{i}.maxcoef(1), rezMax,'UniformOutput',0));
%         lengthPx = cell2mat(cellfun(@(x) x{i}.lengthMatch,rezMax,'UniformOutput',0));
%         pos = cell2mat(cellfun(@(x) x{i}.pos(1), rezMax,'UniformOutput',0));
%         stretch =  cell2mat(cellfun(@(x) x(i),bestBarStretch,'UniformOutput',0));
%         for j=1:length(pos)
%            if pos(j)<= 0
%                pos(j) = pos(j)+thrLen(j);
%            end
%         end
%         maxcc = maxccoef;
%         
%         [d,name,ext] = fileparts(barGen{i}.name);
% 
%         N = matlab.lang.makeValidName(name);
%         
%         T2 = table(maxcc',lengthPx', pos',stretch' ,'VariableNames',{N,strcat(['len_'  num2str(i)]),strcat(['pos_'  num2str(i)]),strcat(['stretch_'  num2str(i)])});
%         T = [T T2];
%     end
%     disp('Saving ccvals table');
%     import CBT.Hca.Export.export_cc;
%     export_cc(T, matDirpath);      
% 
% 
% end
% 
end