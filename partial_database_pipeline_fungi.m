function [] = partial_database_pipeline_fungi()
    % run against database_fungi for semi-realistic barcodes

    thryName = '/export/scratch/albertas/download_dump/fungi/single/theoryOutput/theoryGen_0.34_110_300_0_2024-05-06_12_56_25_session.mat';
    thryName = 'C:\Users\Lenovo\postdoc\Chalmers\8_other\yeast_mappability\theoryGen_0.34_110_300_0_2024-05-06_12_56_25_session.mat';

    import Core.load_theory_structure;
    [theoryStruct,sets] = load_theory_structure(0.2,[],[],thryName);

    for i=1:length(theoryStruct)
        theoryStruct(i).rawBitmask = ones(1,length(  theoryStruct(i).rawBarcode ));
    end
    % todo: maybe don't add extra psf?
    w = 200;

    % remove empty theories
    theoryStruct(arrayfun(@(x) isempty(theoryStruct(x).rawBarcode),1:length(theoryStruct))) = [];
    
    %%

    import Core.Discriminative.extract_species_name;
    [uniqueSpeciesNames,idSpecies] = Core.Discriminative.extract_species_name({theoryStruct.name});
    
    seqNameToTest = 'Candida albicans';
    idSeq = find(cellfun(@(x) isequal(x,seqNameToTest),uniqueSpeciesNames));
    allSeq = find(idSpecies==idSeq);

    totalLength = sum([theoryStruct(allSeq).length]); % how many pixels these contain

    sum([theoryStruct(:).length])

    %% Create a new vector from relevant theories
    thryToCalc = theoryStruct(allSeq);

    % blend forward and reverse
    tempCell = cell(1, 2*numel(thryToCalc)); % 1 is forward, 3 is reverse, 2 and 4 is nan
    tempCell(1:2:end) = {thryToCalc(:).rawBarcode};
%     tempCell(3:4:end) = cellfun(@(x) fliplr(x),{thryToCalc(:).rawBarcode},'un',false);
    tempCell(2:2:end) = {NaN};
    % Concatenate the cell array into a single vector
    vecConcat = cat(2, tempCell{:});

    % convert vec to integers for simpler calculation
    newvecA = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);

    writematrix(newvecA,'barA.txt','Delimiter',' ');

    nanIndices = find(isnan(newvecA));

    % Preallocate memory
    indexes = cell(1,numel(nanIndices)/2);
    % Split the vector based on NaN delimiters
    for i = 1:numel(nanIndices) % every two, since doesn't matter if we look at forward or reverse
        if i == 1
            indexes{i} = [1 nanIndices(i)-1];
        else
            indexes{i}  = [nanIndices(i-1)+1 nanIndices(i)-1];
        end
    end


    %% For rest of the theories
    wmin = 600;

    thryToCalc = theoryStruct(idSpecies~=idSeq);
    thryToCalc = thryToCalc(find([thryToCalc.length]>wmin));

    tempCell = cell(1, 4*numel(thryToCalc)); % 1 is forward, 3 is reverse, 2 and 4 is nan
    tempCell(1:4:end) = {thryToCalc(:).rawBarcode};
    tempCell(3:4:end) = cellfun(@(x) fliplr(x),{thryToCalc(:).rawBarcode},'un',false);
    tempCell(2:2:end) = {NaN};
    % Concatenate the cell array into a single vector
    vecConcat = cat(2, tempCell{:});

    % convert vec to integers for simpler calculation
    newvecB = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);

    writematrix(newvecB,'barB.txt','Delimiter',' ');

    %%
    

%     strjoin([theoryStruct(1:2).rawBarcode])
        numWorkers = 30;
       com= strcat(['SCAMP --window=' num2str(wmin) ' --input_a_file_name='...
           fullfile(pwd,'barA.txt') ' --input_b_file_name=' ...
           fullfile(pwd,'barB.txt') ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
           ' --output_a_file_name=' 'bar_mp' ...
           ' --output_a_index_file_name=' 'bar_index']);

    tic
    [a,val ] = system(com);
    toc


    %%
    
    fid = fopen('bar_mp');
    raw2 = textscan(fid, '%s ');
    fclose(fid);
    nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
    mp1 = nan(length(nonanValues),1);
    mp1(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');

    for i=1:length(indexes)
        mp1(indexes{i}(1):indexes{i}(2))
    end

    figure,tiledlayout(3,2)
    nexttile([1 2])
    imagesc([mp1>0.85]');colormap(gray)
    cmap = [.7 .7 .7 %// light gray
        .5 .5 .5] %// white
    colormap(cmap)
    axis off
    nexttile
    xticks([0 5 10])

    % colorbar('Ytick',[.25 .75],'Yticklabel',[0 1]) %// only two values in colorbar

%%
% % 
% % lenThr = 10000;
% % figure,tiledlayout(3,1)
% % nexttile
% % plot(newvec(1:lenThr))
% % title('Theoretical barcodes')
% % nexttile
% % plot(mp1(1:lenThr))
% % title('Local score, w=200 (~400kb)')
% % nexttile
% % plot(mpI(1:lenThr))
% % title('Local index')
% 
% 
% nanIndices = nanIndices(nanIndices<length(mp1)-wmin);
% 
% mapMatrix = nan(numel(nanIndices)/2,numel(nanIndices)/2);
% for i=1:numel(nanIndices)/2
%     i
%     % curM =  mp1(indexes{i}(1):min(end,indexes{i}(2)));
%     curElements = mp1(indexes{i}(1):min(end,indexes{i}(2)));
%     curElements = curElements(curElements~=-1);
%     barcodesMappedTo = barIndex(curElements+1);
%     [vals,rep] = unique(barcodesMappedTo);
%     mapMatrix(i,vals) = rep/length(barcodesMappedTo);
% end
% 
% 
% 
% %%
% 
% % 
% % 
% %     import Core.rescale_barcode_data;
% %     bG = [];
% %     bG{1} = theoryStruct(allSeq(1));
% % %     bG{1}.rawBitmask = theoryStruct(allSeq(1)).rawBitmask(1:wmin);
% % % 
% % %     bG{1}.rawBarcode = theoryStruct(allSeq(1)).rawBarcode(1:wmin);
% % %     bG{1}.rawBitmask = theoryStruct(allSeq(1)).rawBitmask(1:wmin);
% % 
% %     [barGen] = rescale_barcode_data(bG,1);
% % 
% % 
% % 
% % import CBT.Hca.Core.Comparison.hca_compare_distance;
% % import Core.export_coefs_local;
% % 
% % %     display(['Running w = ',num2str(sets.w)]);
% %     % only local lengths for which all length re-scaled versions passes the
% %     % threshold
% % %     passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),barGen)*sets.theory.stretchFactors(end) >= sets.w);
% % 
% % sets.w = 600;
% % sets.comparisonMethod = 'mpnan'; % 'mpnan'
% % 
% % % Compare to sequences which are not seqNameToTest
% % 
% % [rezMaxMP] = hca_compare_distance(barGen, theoryStruct, sets );
% 
% 
% 
%     %% Database self-similarity. A bit slow if to run all pairwise. So want to make the theory smaller
%     % comparisons. If we run long vs long, and keep only barcodes with
% %     % subbarcodes score < thresh, maybe more accurate?
% %     tic % test:
% %     [oS] = calc_overlap_mp(theoryStruct(1:20),1, w, 'test');
% %     toc
% 
% 
%     %% Create a new vector from both theories and their reverse complements.
%     thryToCalc = theoryStruct(1:200);
% 
%     % vec = [thryToCalc(:).rawBarcode nan];
%     tempCell = cell(1, 4*numel(thryToCalc)); % 1 is forward, 3 is reverse, 2 and 4 is nan
%     tempCell(1:4:end) = {thryToCalc(:).rawBarcode};
%     tempCell(3:4:end) = cellfun(@(x) fliplr(x),{thryToCalc(:).rawBarcode},'un',false);
%     tempCell(2:2:end) = {NaN};
%     % Concatenate the cell array into a single vector
%     vecConcat = cat(2, tempCell{:});
% 
%     % convert vec to integers for simpler calculation
%     newvec = round((vecConcat-min(vecConcat))/(max(vecConcat)-min(vecConcat))*256);
% 
%     %%
%     nanIndices = find(isnan(newvec));
% 
%     % Preallocate memory
%     barIndex = zeros(1,length(newvec));
%     indexes = cell(1,numel(nanIndices)/2);
%     % Split the vector based on NaN delimiters
%     for i = 1:2:numel(nanIndices) % every two, since doesn't matter if we look at forward or reverse
%         if i == 1
%             barIndex(1:nanIndices(i+1)-1) = i*ones(1, nanIndices(i+1)-1);
%             indexes{i} = [1 nanIndices(i+1)-1];
%         else
%             barIndex(nanIndices(i-1)+1:nanIndices(i+1)-1) = ceil(i/2)*ones(1, nanIndices(i+1)-1-(nanIndices(i-1)+1)+1);
%             indexes{ceil(i/2)}  = [nanIndices(i-1)+1 nanIndices(i+1)-1];
%         end
%     end
%     
%     % Don't need the last element anymore
%     % % Add the last part of the vector if NaN is not the last element
%     % if nanIndices(end) ~= numel(barIndex)
%     %     barIndex(nanIndices(end)+1:end) = i*ones(1,numel(barIndex)-nanIndices(end));
%     %     indexes{numel(nanIndices)+1} = [nanIndices(end)+1 numel(barIndex)];
%     % end
% 
% %%
% 
%     %
% %     newvec(isnan(vecConcat)) = nan;
%     writematrix(newvec,'bar.txt','Delimiter',' ')
% %     strjoin([theoryStruct(1:2).rawBarcode])
% numWorkers = 30;
%        com= strcat(['SCAMP --window=' num2str(w) ' --input_a_file_name='...
%            fullfile(pwd,'bar.txt') ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
%            ' --output_a_file_name=' 'bar_mp' ...
%            ' --output_a_index_file_name=' 'bar_index']);
% 
%     tic
%     [a,val ] = system(com);
%     toc
% 
% % % comp time
% % x = [100 200 400 1000];
% % y = [1.5 5 16.3 100];
% % coefs = polyfit(x,y,2);
% % figure,plot(x,y);hold on
% % xx=100:1000;
% % timefun = @(xx) coefs(1)*xx.^2+coefs(2)*xx+coefs(3);
% % plot(xx,timefun(xx))
% 
% %
% 
% fid = fopen('bar_mp');
% raw2 = textscan(fid, '%s ');
% fclose(fid);
% nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
% mp1 = nan(length(nonanValues),1);
% mp1(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%    
% mpI = importdata('bar_index');
% 
% lenThr = 10000;
% figure,tiledlayout(3,1)
% nexttile
% plot(newvec(1:lenThr))
% title('Theoretical barcodes')
% nexttile
% plot(mp1(1:lenThr))
% title('Local score, w=200 (~400kb)')
% nexttile
% plot(mpI(1:lenThr))
% title('Local index')
% 
% 
% nanIndices = nanIndices(nanIndices<length(mp1)-w);
% 
% %% todo: convert to similarity matrix, using the mpI. Create vector with bar index
% mapMatrix = nan(numel(nanIndices)/2,numel(nanIndices)/2);
% for i=1:numel(nanIndices)/2
%     i
%     % curM =  mp1(indexes{i}(1):min(end,indexes{i}(2)));
%     curElements = mpI(indexes{i}(1):min(end,indexes{i}(2)));
%     curElements = curElements(curElements~=-1);
%     barcodesMappedTo = barIndex(curElements+1);
%     [vals,rep] = unique(barcodesMappedTo);
%     mapMatrix(i,vals) = rep/length(barcodesMappedTo);
% end
% 
% 
% pyoIds = find(idSpecies==3284);
% 
% figure,imagesc(mapMatrix(pyoIds,pyoIds));colormap(gray);colorbar
% title('S.Pyogenes local similarity matrix')
% 
% ecoliId = find(cellfun(@(x) ~isempty(strfind('Escherichia coli',x)),uniqueSpeciesNames));
% ecoliIds = find(idSpecies==ecoliId);
% 
% figure,imagesc(mapMatrix(pyoIds,ecoliIds));colormap(gray);colorbar
% title('S.Pyogenes vs E.Coli local similarity matrix')
% 
% 
% %% now check for which vectors has PCC larger than pccThresh. We can also find corresponding theories that they match
% % best to. Then one of them should be removed and one kept.
% 
% pccThresh = 0.9;
% outVec = mp1 < pccThresh;
% % Preallocate cell array to store splitted vectors
% splittedVectors = zeros(1,length(indexes));
% 
% % Split the vector based on NaN delimiters
% for i = 1:length(indexes)
%     splittedVectors(i) =  mean(outVec(indexes{i}(1):min(end,indexes{i}(2))));
% end
%     % if i == 1
%     %     splittedVectors(i) = mean(outVec(1:nanIndices(i)-1));
%     % else
%     %     splittedVectors(i) = mean(outVec(nanIndices(i-1)+1:nanIndices(i)-1));
%     % end
% % end
% 
% % 
% % % Add the last part of the vector if NaN is not the last element
% % if nanIndices(end) ~= numel(outVec)
% %     splittedVectors(end+1) = mean(outVec(nanIndices(end)+1:end));
% % end
% 
% % Repetitive sequences: challenge is which sequence to keep. Maybe have to
% % look into mpI
% 
% figure,histogram(splittedVectors)
% title(['Theories locally matching somewhere with PCC <', num2str(pccThresh)])
% xlabel(['Percentage of theories with PCC < ',num2str(pccThresh)])
% 
% 
% 
% figure,histogram(splittedVectors(pyoIds),'Normalization','pdf');
% hold on
% histogram(splittedVectors(ecoliIds),'Normalization','pdf');
% title(['Theories locally matching somewhere with PCC < ', num2str(pccThresh)])
% xlabel(['Percentage of theories with PCC < ', num2str(pccThresh)])
% legend({'S.Pyogenes','E.Coli'})
% 
% 
% %%
% thrid1 = 3538;
% thrid2 = 3284;
% 
% id1 = find(idSpecies==thrid1);
% id2 = find(idSpecies==thrid2);
% 
% figure,imagesc(mapMatrix(id1,id2));colormap(gray);colorbar
% title('Local similarity matrix')
% 
%        %%
% 
%     timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
%     foldSynth=strcat('output',timestamp);
% 
%     Ntheories = 10;
%     [~,~] =mkdir(foldSynth);
% 
% 
%     % have to save separately..
%     [namesBar, stridx,barStruct] = Core.save_bars_rescaled_txt(barGen,sF,foldSynth);
%     % save one long theory
%     [names2, baridx2] = Core.save_long_theory(theoryStruct(1:Ntheories),'barcoli');
% 
% tS{1}.filename = names2{1};
% 
% numWorkers = 30;
% tic
% sF = 1;
% sets.w = 300;
% import Core.compare_mp_multi_theories_fast
% % for i=1:length(barGen)
% %     [compStrOut] = compare_mp_multi_theories(barGen(i),tS,sets.theory.stretchFactors,sets.w,numWorkers);
% % end
% % toc
% 
% barGen = theoryStruct(1);
% barGen(1).rawBarcode = barGen(1).rawBarcode(1:300);
% barGen(1).rawBitmask = logical(ones(1,length(barGen(1).rawBarcode)));
% 
% tic
% [overlapStruct2] = compare_mp_multi_theories_fast(barGen, theoryStruct, 0.95:0.025:1.05, w, numWorkers);
% toc
% 
% % PLOT. First calculate full overlap
% ix = 90;
% theoryStruct(ix).rawBitmask = logical(ones(1,length(theoryStruct(ix).rawBarcode ) ));
% import Core.get_full_overlap_score;
% [overlapStruct2(1,ix).fullscore,overlapStruct2(1,ix).overlaplen, overlapStruct2(1,ix).lenB , overlapStruct2(1,ix).lenA,overlapStruct2(1,ix).partialScore,...
%      overlapStruct2(1,ix).partialLength] = get_full_overlap_score(overlapStruct2(1,ix).pA,overlapStruct2(1,ix).pB,...
%      overlapStruct2(1,ix).bestBarStretch, overlapStruct2(1,ix).or,[barGen(1) theoryStruct(ix)],w);
% 
% 
% %f= figure;
% % tiledlayout(2,1)
% import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
% pair_evaluation_with_ground_truth_plot([barGen(1) theoryStruct(ix)], [overlapStruct2(1,1) overlapStruct2(1,ix)],1,2,[],10000);
% 
% 
% 
% 
% end
% 
end