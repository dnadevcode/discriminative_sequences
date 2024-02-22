function [] = partial_database_pipeline()
    % run against partial database




    % loading theory to specific nm/bp
    tic
    import Core.load_theory_structure;
    thryFileIdx = 1; % todo: only pass single
    [theoryStruct,sets] = load_theory_structure(0.05,thryFileIdx);
    toc
    for i=1:length(theoryStruct)
        theoryStruct(i).rawBitmask = ones(1,length(  theoryStruct(i).rawBarcode ));
    end
    % todo: maybe don't add extra psf?

    %% Database self-similarity. A bit slow if to run all pairwise. So want to make the theory smaller
    % comparisons. If we run long vs long, and keep only barcodes with
    % subbarcodes score < thresh, maybe more accurate?
    tic
    [oS] = calc_overlap_mp(theoryStruct(1:50),1, 200, 'test');
    toc

    vec = [theoryStruct(:).rawBarcode nan];
    tempCell = cell(1, 2*numel(theoryStruct)-1);
    tempCell(1:2:end) = {theoryStruct(:).rawBarcode};
    tempCell(2:2:end) = {NaN};
    % Concatenate the cell array into a single vector
    vecConcat = cat(2, tempCell{:});


    % convert vec to integers
    newvec = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);
%     newvec(isnan(vecConcat)) = nan;
    writematrix(newvec,'bar.txt','Delimiter',' ')
%     strjoin([theoryStruct(1:2).rawBarcode])
w = 200;
numWorkers = 32;
       com= strcat(['SCAMP --window=' num2str(w) ' --input_a_file_name='...
           '/export/scratch/albertas/data_temp/bargrouping/PARTIAL_DB_DATA/bar.txt' ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=' 'bar_mp' ...
           ' --output_a_index_file_name=' 'bar_index']);

tic
[a,val ] = system(com);
toc

fid = fopen('bar_mp');
raw2 = textscan(fid, '%s ');
fclose(fid);
nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
mp1 = nan(length(nonanValues),1);
mp1(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
   
figure,plot(mp1)


nanIndices = find(isnan(newvec));
nanIndices = nanIndices(nanIndices<length(mp1)-w);

pccThresh = 0.9;
outVec = mp1 < pccThresh;
% Preallocate cell array to store splitted vectors
splittedVectors = zeros(1, numel(nanIndices));

% Split the vector based on NaN delimiters
for i = 1:numel(nanIndices)
    if i == 1
        splittedVectors(i) = mean(outVec(1:nanIndices(i)-1));
    else
        splittedVectors(i) = mean(outVec(nanIndices(i-1)+1:nanIndices(i)-1));
    end
end

% Add the last part of the vector if NaN is not the last element
if nanIndices(end) ~= numel(outVec)
    splittedVectors(end+1) = mean(outVec(nanIndices(end)+1:end));
end

% Repetitive sequences: challenge is which sequence to keep. Maybe have to
% look into mpI



       %%

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    foldSynth=strcat('output',timestamp);

    Ntheories = 10;
    [~,~] =mkdir(foldSynth);


    % have to save separately..
    [namesBar, stridx,barStruct] = Core.save_bars_rescaled_txt(barGen,sF,foldSynth);
    % save one long theory
    [names2, baridx2] = Core.save_long_theory(theoryStruct(1:Ntheories),'barcoli');

tS{1}.filename = names2{1};

numWorkers = 30;
tic
sF = 1;
sets.w = 300;
import Core.compare_mp_multi_theories_fast
% for i=1:length(barGen)
%     [compStrOut] = compare_mp_multi_theories(barGen(i),tS,sets.theory.stretchFactors,sets.w,numWorkers);
% end
% toc

barGen = theoryStruct(1);
barGen(1).rawBarcode = barGen(1).rawBarcode(1:300);
barGen(1).rawBitmask = logical(ones(1,length(barGen(1).rawBarcode)));

tic
[overlapStruct2] = compare_mp_multi_theories_fast(barGen, theoryStruct, 0.95:0.025:1.05, w, numWorkers);
toc

% PLOT. First calculate full overlap
ix = 90;
theoryStruct(ix).rawBitmask = logical(ones(1,length(theoryStruct(ix).rawBarcode ) ));
import Core.get_full_overlap_score;
[overlapStruct2(1,ix).fullscore,overlapStruct2(1,ix).overlaplen, overlapStruct2(1,ix).lenB , overlapStruct2(1,ix).lenA,overlapStruct2(1,ix).partialScore,...
     overlapStruct2(1,ix).partialLength] = get_full_overlap_score(overlapStruct2(1,ix).pA,overlapStruct2(1,ix).pB,...
     overlapStruct2(1,ix).bestBarStretch, overlapStruct2(1,ix).or,[barGen(1) theoryStruct(ix)],w);


%f= figure;
% tiledlayout(2,1)
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_plot([barGen(1) theoryStruct(ix)], [overlapStruct2(1,1) overlapStruct2(1,ix)],1,2,[],10000);




end

