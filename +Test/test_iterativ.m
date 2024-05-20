% iterative re-scaling
% have initial consensus

cInitial = cGenAll{ix};

ccthresh = 0;
NN = 20;
consensusNew = cell(1,NN+1);
consensusNew{1} = consensus;

for jj=1:5
    jj

% with re-scaling
barC = barcodeGen(cInitial.idx);
initialStretch = cellfun(@(x) x.bestBarStretch,cInitial.comparisonStruct);

thr = [];
thr(1).rawBarcode = consensusNew{jj}(~isnan(consensusNew{jj}));
thr(1).rawBitmask = ~isnan(thr(1).rawBarcode) ;
thr(1).isLinearTF = 1;

sF = 0.98:0.01:1.02;

import Core.rescale_barcode_data;
[barCRescaled] = rescale_barcode_data(barC,sF,initialStretch);

sets.comparisonMethod = 'mass_pcc';
sets.w = nan;
import CBT.Hca.Core.Comparison.hca_compare_distance;
[rezMaxMP] = hca_compare_distance(barCRescaled, thr, sets );
    
% create overlap table
%     info: {'1)PCC, 2) Pos, 3) Or, 4) SecondPos 5) Len '}
vals = double(zeros(size(rezMaxMP{1}{1},1),4));
score = zeros(size(rezMaxMP{1}{1},1),1);
for i=1:length(barCRescaled)
    [maxV,maxPos] = max(rezMaxMP{1}{1}(i,:));
    posStart = double(rezMaxMP{1}{2}(i,maxPos));
    posStop = double(rezMaxMP{1}{2}(i,maxPos)+length(barCRescaled{i}.rescaled{maxPos}.rawBarcode)-1);
    posOr = double(rezMaxMP{1}{3}(i,maxPos));
    curSF =  initialStretch(i)*sF(maxPos);
    score(i) = maxV;

%     if (score(i) < ccthresh) && jj>1 || curSF<sF(1) || curSF > sF(end)
%         vals(i,:) = valsprev(i,:);% [cInitial.comparisonStruct{i}.pos cInitial.comparisonStruct{i}.pos+cInitial.comparisonStruct{i}.lengthMatch-1 cInitial.comparisonStruct{i}.or cInitial.comparisonStruct{i}.bestBarStretch]; % all these positions are w.r.t. the consensus
%     else
        vals(i,:) = [posStart posStop posOr curSF]; % all these positions are w.r.t. the consensus
%     end

%     stats.posdif(i) = cInitial.comparisonStruct{i}.pos-posStart;
%     stats.ordif(i) = isequal(cInitial.comparisonStruct{i}.or,posOr);
%     stats.sfDif(i) = sF(maxPos);
end

vals = vals(score>ccthresh,:);
ids = barIslands{ix}(score>ccthresh);

import Plot.islandsPosStruct;
cGenNew.comparisonStruct = islandsPosStruct({vals},{ids});
cGenNew.idx = cInitial.idx(score>ccthresh);
import Core.barcode_island_consensus;
multiDimBarNew = barcode_island_consensus(barcodeGen,{cGenNew}, 1, 300);

consensusNew{jj+1} = nanmean(multiDimBarNew);

valsprev = vals;
cInitial = cGenNew;
stats.score{jj} = score;

end

%%  consensus vs theory
% First generate theory:
w = 300;% minimum overlap length
% sF = 0.9:0.025:1.1;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

fastas = {'018_final_polish.fasta','DA32087.fasta','DA68335.fasta'};
theoryIdxs = {1, 3, 2, nan, 2 , 2, 3}; % known theory indexes
dataSetIdx = 3; % 5,... test datasets (part of the data)

nmPerPx = 110;
nmbp = 0.25; % ? 

% quickly calculate theory for given single fasta file
[theoryStructRev,theoryStruct,bT] = prep_thry_for_local_comp(fastas(theoryIdxs{dataSetIdx}), nmbp, nmPerPx, 1);

%%
sF = 0.9:0.01:1.1;
theoryStruct.rawBitmask = zeros(1,length(theoryStruct.rawBarcode));
import Core.rescale_barcode_data;
[tsRescaled] = rescale_barcode_data({theoryStruct},sF);

import SignalRegistration.masked_multidim_pcc_corr;
maxCoefsConsensus = zeros(1,length(consensusNew));  
for k=1:length(consensusNew)
    xcorrsBest = cell(1,length(tsRescaled{1}.rescaled));
    for i=1:length(tsRescaled{1}.rescaled)
        [ xcorrsBest{i}, ~ ] = masked_multidim_pcc_corr( consensusNew{k}',tsRescaled{1}.rescaled{i}.rawBarcode',~isnan(consensusNew{k})',~isnan(tsRescaled{1}.rescaled{i}.rawBarcode)',300);
    end
    
    zMethod = 'meanpcc';
    scoreSingleConsensus = cell(1,length(xcorrsBest));
    for i=1:length(xcorrsBest)
        [maxScore,maxPos] = max(xcorrsBest{i},[],1);
        % first take max over 3rd dim
        [scoreMaxOr,scoreMaxOrPos] = max(xcorrsBest{i},[],3);
        switch zMethod
            %             case 'stouffer'
            %             case 'fisher'
            case 'meanpcc'
                scoreSingleConsensus{i} = mean(scoreMaxOr,2);
            otherwise
                error('No such method for combining p-values')
        end
    end
    
    maxCoefsConsensus(k) = max(cellfun(@(x) max(x),scoreSingleConsensus));

end

%
figure,plot(maxCoefsConsensus);
xlabel('Iteration')
ylabel('Max PCC consensus vs theory')