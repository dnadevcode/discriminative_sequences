function [] = partial_database_pipeline()
    % run against partial database





    tic
    import Core.load_theory_structure;
    thryFileIdx = 1; % todo: only pass single
    [theoryStruct,sets] = load_theory_structure(0.1,thryFileIdx);
    toc
    % todo: maybe don't add extra psf?

    % Database self-similarity


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

