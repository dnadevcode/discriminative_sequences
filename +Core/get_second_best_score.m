function [] = get_second_best_score(idSpecies,idSeq,theoryStruct,wmin,numWorkers)

%% For rest of the theories
restOfTheories = find(idSpecies~=idSeq);
thryToCalc = theoryStruct(idSpecies~=idSeq);

% thryToCalc = theoryStruct(1:100);
tempCell = cell(1, 4*numel(thryToCalc)); % 1 is forward, 3 is reverse, 2 and 4 is nan
tempCell(1:4:end) = {thryToCalc(:).rawBarcode};
tempCell(3:4:end) = cellfun(@(x) fliplr(x),{thryToCalc(:).rawBarcode},'un',false);
tempCell(2:2:end) = {NaN};
% Concatenate the cell array into a single vector
vecConcat = cat(2, tempCell{:});

% convert vec to integers for simpler calculation
newvecB = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);

nanIndices = find(isnan(newvecB));

% Preallocate memory
indexesT = cell(1,numel(nanIndices)/2);
% Split the vector based on NaN delimiters
for i = 1:numel(nanIndices) % every two, since doesn't matter if we look at forward or reverse
    if i == 1
        indexesT{i} = [1 nanIndices(i)-1];
    else
        indexesT{i}  = [nanIndices(i-1)+1 nanIndices(i)-1];
    end
end


writematrix(newvecB','barB.txt','Delimiter',' '); 

% numWorkers = 4;

com= strcat(['SCAMP --window=' num2str(wmin) ' --input_a_file_name='...
    fullfile(pwd,'barA.txt') ' --input_b_file_name=' ...
    fullfile(pwd,'barB.txt') ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
    ' --output_a_file_name=' 'bar_mp' ...
    ' --output_a_index_file_name=' 'bar_index']);

tic
[a,val ] = system(com);
toc
end

