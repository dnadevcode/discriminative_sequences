function [] = disc_run_simulation(thryName,nmbp)

% add to path if not added
addpath(genpath('/home/avesta/albertas/reps/hca'))
addpath(genpath('/home/avesta/albertas/reps/discriminative_sequences'))


if nargin < 1
    % run against database_fungi for semi-realistic barcodes
    %     thryName = '/export/scratch/albertas/download_dump/fungi/single/theoryOutput/theoryGen_0.34_110_300_0_2024-05-06_12_56_25_session.mat';
    thryName = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\yeast_mappability\theoryGen_0.34_110_300_0_2024-05-06_12_56_25_session.mat';

end

if nargin < 2 || isempty(nmbp)
    nmbp = 0.22;
end


nmpx = 110; %nm/px
psf = 300; %nm



import Core.load_theory_structure;
[theoryStruct,sets] = load_theory_structure(nmbp,[],[],thryName);

for i=1:length(theoryStruct)
    theoryStruct(i).rawBitmask = ones(1,length(  theoryStruct(i).rawBarcode ));
end
% todo: maybe don't add extra psf?
w = 200;
wmin = 600;


% remove empty theories and the ones shorter than some minimum length wmin
theoryStruct(arrayfun(@(x) isempty(theoryStruct(x).rawBarcode),1:length(theoryStruct))) = [];
theoryStruct = theoryStruct(find([theoryStruct.length]>wmin));


%% 

import Core.Discriminative.extract_species_name;
[uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name({theoryStruct.name});

% length(find(cellfun(@(x) contains(x,'Candida'),{theoryStruct.name})))/length({theoryStruct.name})
seqNameToTest = 'Candida albicans';
%seqNameToTest = '[Candida] auris';

idSeq = find(cellfun(@(x) isequal(x,seqNameToTest),uniqueSpeciesNames));
allSeqA = find(idSpecies==idSeq);
numel(allSeqA)
{theoryStruct(allSeqA).name}' % chromosomes are ordered like this


seqNames = cellfun(@(u) u(1),cellfun(@(z) strsplit(z,'chromosome'), cellfun(@(x) x(2),arrayfun(@(y) regexp(theoryStruct(allSeqA(y)).name,' ','split','once'),1:length(allSeq),'UniformOutput',false)),'un',false))
[uniqueSeqNames, ~, idSeq] = unique( seqNames ) ; % todo: if some species are considered the same, include that here

chrid = 1;
allSeq  =find(idSeq==chrid);


% chr = cellfun

CCmax = 0.7;

%% Create a new vector from relevant theories
thryToCalc = theoryStruct(allSeq);

% thryToCalc = theoryStruct(allSeq(1:7));

% only forward 
tempCell = cell(1, 2*numel(thryToCalc)); % 1 is forward 2 is nan
tempCell(1:2:end) = {thryToCalc(:).rawBarcode};
tempCell(2:2:end) = {NaN};

% Concatenate the cell array into a single vector
vecConcat = cat(2, tempCell{:});
% convert vec to integers for simpler calculation
newvecA = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);
writematrix(newvecA,'barC.txt','Delimiter',' ');


% simple rand procedure: first fix the sRand based on similarity to the
% beginning of the sequence. We need to do this for a few random regions
% along the genome in order to get close to average CCmax in practice
N = 100;
thrId = 1;
sRandV = zeros(1,N);
for i=1:N
    rb = imgaussfilt(normrnd(0,1,1,w),psf/nmpx);
    stPos = randi(length(tempCell{thrId})-w+1,1); % starting pos on thry;
    % stPos = 5000;
    x = 0.1;
    fz = @(x) zscore( tempCell{thrId}(stPos:stPos+w-1),1)*zscore( tempCell{thrId}(stPos:stPos+w-1)+x*rb,1)'/w- CCmax;
    sRandV(i) = fzero(@(x) fz(x), [0 1]);
end
sRand = mean(sRandV,'omitnan'); % sRand for this level of noise

% now add to tempCell
for i=1:2:numel(tempCell)
    tempCell{i} = tempCell{i}+sRand*imgaussfilt(normrnd(0,1,1,length(tempCell{i})),psf/nmpx);
end


% Concatenate the cell array into a single vector
vecConcat = cat(2, tempCell{:});

% convert vec to integers for simpler calculation
newvecA = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);

writematrix(newvecA,'barA.txt','Delimiter',' ');

nanIndices = find(isnan(newvecA));

% Preallocate memory
indexes = cell(1,numel(nanIndices));
% Split the vector based on NaN delimiters
for i = 1:numel(nanIndices) % every two, since doesn't matter if we look at forward or reverse
    if i == 1
        indexes{i} = [1 nanIndices(i)-1];
    else
        indexes{i}  = [nanIndices(i-1)+1 nanIndices(i)-1];
    end
end

figure,plot(tempCell{1}(1:600))


numWorkers = 30;
import Core.get_second_best_score;
[indexesT,restOfTheories] = get_second_best_score(idSpecies,idSeq,theoryStruct,wmin,numWorkers);
%


com= strcat(['SCAMP --window=' num2str(wmin) ' --input_a_file_name='...
    fullfile(pwd,'barA.txt') ' --input_b_file_name=' ...
    fullfile(pwd,'barC.txt') ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
    ' --output_a_file_name=' 'bar_mpC' ...
    ' --output_a_index_file_name=' 'bar_indexC']);

tic
[a,val ] = system(com);
toc

end

