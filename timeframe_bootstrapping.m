function outputArg1,outputArg2] = timeframe_bootstrapping(barGenRun)
% bootstrapping using aligned kymo time-frames

alignedKymo = barGenRun{1}.alignedKymo; % should not have nans in the middle.. check the generation
leftEdgeIdxs =  barGenRun{1}.leftEdgeIdxs;
rightEdgeIdxs =  barGenRun{1}.rightEdgeIdxs;
rawBitmask = barGenRun{1}.rawBitmask;
numWindows = size(alignedKymo,1);

import DBM4.gen_barcode_data;
%     barcodeGen =  gen_barcodes_from_kymo( kymoStructs{idFold}, sets,sets.maxLen);
 

if ~exist('s','var')
    s = RandStream('mlfg6331_64'); 
end

NN = 1000;
%     pccScore = zeros(1,NN);
barcodeGenBootstrapped = cell(1,NN);
stdmat = zeros(NN,length(barGenRun{1}.lE:barGenRun{1}.rE));
meanmat = zeros(NN,length(barGenRun{1}.lE:barGenRun{1}.rE));

for i=1:NN
    y = datasample(s,1:numWindows,numWindows,'Replace',true);
    kymoBootstrapped = alignedKymo(y,:);
    [barcodeGenBootstrapped{i}] = gen_barcode_data(kymoBootstrapped,leftEdgeIdxs(y), rightEdgeIdxs(y),1);
    barcodeGenBootstrapped{i}.rawBarcode =   barcodeGenBootstrapped{i}.rawBarcode(barGenRun{1}.lE:barGenRun{1}.rE); % check 
    barcodeGenBootstrapped{i}.rawBitmask = rawBitmask;
    stdmat(i,:) = nanstd(kymoBootstrapped(:,barGenRun{1}.lE:barGenRun{1}.rE),1);
    meanmat(i,:) = nanmean(kymoBootstrapped(:,barGenRun{1}.lE:barGenRun{1}.rE),1);

end
    % 

    f = figure
    tiledlayout(2,1)
    nexttile
    hold on
    cellfun(@(x) plot(zscore(x.rawBarcode),'color',   [0.5 0.8 0.8],'LineWidth',0.1),barcodeGenBootstrapped)
    nexttile
    super_quick_plot(1,barGenRun,rezOut{1},theoryStruct,f)


barGenRun2 =barGenRun;
barGenRun2{1}.rawBarcode = barGenRun2{1}.alignedKymo(1,barGenRun{1}.lE:barGenRun{1}.rE);
barGenRun2{1}.rawBarcode(isnan(barGenRun2{1}.rawBarcode)) = nanmean(barGenRun2{1}.rawBarcode);
        super_quick_plot(1,barGenRun2,rezOut{1},theoryStruct);

end

