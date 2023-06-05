% theory barcodes location
% theory_names = dir('/proj/snic2022-5-384/users/x_albdv/data/chromseq/*.fasta');

% N = 100; % number theory barcodes to take

data = theory_names(1:100);%(1:N);

% settings file
setsFile= 'theory_settings_parallel.txt';
import CBT.Hca.Import.import_settings;
[sets] = import_settings(setsFile);

sets.addgc = 0;
sets.resultsDir = fullfile(pwd,sets.resultsDir);

% fasta files
dataF = arrayfun(@(x) fullfile(x.folder,x.name),data,'UniformOutput',false);

% sets.theoryGen

fd = fopen('theories_parallel.txt','w');
for i=1:length(data)
fprintf(fd,'%s\n',dataF{i});
end
fclose(fd);

sets.fastas = 'theories_parallel.txt';
sets.skipBarcodeGenSettings=1;
sets.skipChangeBpNmRatio=1;
sets.theoryGen.pixelWidth_nm = 110;
mkdir(sets.resultsDir);

% 2 ways to calculate. om_theory todo: use this (since allows to include
% settings
% tic
% [t,matFilepathShort,theoryStruct, sets,theoryGen] = HCA_om_theory_parallel(1,0.25,sets);
% toc

% directly calculate
tic
[theoryStructRev,theoryStruct,barcodeGen] = prep_thry_for_local_comp(dataF, 0.34, sets.theoryGen.pixelWidth_nm, 1);
toc

tic
minLen = 250;
MP = cell(1,length(minLen));
MPI = cell(1,length(minLen));
for ii=1:length(minLen)
    [MP{ii}, MPI{ii}] = local_compare_thry(minLen(ii),[], [], [], [], 30,theoryStructRev(1));
end
toc

mpMax = cellfun(@(x) max(x{1}(1:end/2)),MP);

% maybe use prev info of which barcode match to which

% Idea: will need to bitmask regions which have high pcc

% mpMaxEcoli = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,psffac,numWorkers,minLen);

% %% then run MP self-similarity
% % sum(cellfun(@(x) x.length,theoryStruct))
% 
% barStruct = cell2struct([cellfun(@(x) x,theoryGen.theoryBarcodes,'un',false);...
%     cellfun(@(x) x<200,theoryGen.theoryBitmasks,'un',false)]',{'rawBarcode','rawBitmask'},2);
% 
% % todo: include circular
% foldSynth = 'datasingle';
% % save as single txt.
% [names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,[inf],foldSynth);
% 
% %%
% tic
% numWorkers = 30;
% mpI1  = [];
% mp1 = [];
% MIN_OVERLAP_PIXELS = 250;
% %     numWorkers = 4;
% %     for k=1:length(names2)
% k = 1;
% 
%     if ispc
%         command = strcat(['SCAMP --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
%         names2{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
%         ' --output_a_file_name=.\' fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
%         ' --output_a_index_file_name=.\' fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
%     else
%     command = strcat(['~/SCAMP/build/SCAMP --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
%     names2{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
%     ' --output_a_file_name=\' fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
%     ' --output_a_index_file_name=\' fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
%     end
%         
%         [~,message] = system(command);
% toc        
% %% load
% 
% 
%         mpI1{k} = importdata(fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_index' ])));
%     
%     %also get PCC scores.
% %     tic
%     % a bit slower import to make sure imports correctly on win
%     fid = fopen(fullfile(foldSynth,strcat([num2str(k) '_' num2str(k) '_mp' ])));
%     raw2 = textscan(fid, '%s ');
%     fclose(fid);
%     nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
%     mp1{k} = nan(length(nonanValues),1);
%     mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
% %     toc
% 
% 
% %     end
%     %%
%     
%     % maximum values,, 2 max, since two locations should be
%     % similar
%     
%     sortedMP = cellfun(@(x) sort(x,'desc','MissingPlacement','last'),mp1,'un',false);
%     
%     pccValsMax = cellfun(@(x) mean(x(1:2)), sortedMP);
%     
%     valsPassingThresh = (pccValsMax<thresh);
%     
%     
% % end
% 
% % 
% % 
% % function [mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE,MIN_OVERLAP_PIXELS,numWorkers,names2,namesBar)
% % % parpool(c)
% % NN=length(namesBar);
% % import Zeromodel.beta_ev_cdf;
% % 
% % out=strcat('output',timestamp);
% % % mkdir(out);
% % % parpool('local',28)
% % mp1 = cell(1,NN);
% % % mp2 = cell(1,NN);
% % % PKS1=cell(1,NN);
% % % LOCS1=cell(1,NN);
% % % pksUnique1 = cell(1,NN);
% % % pksUniquePos1=cell(1,NN);
% % % barcodePair1=cell(1,NN);
% % % rescaleFactorPair1=cell(1,NN);
% % % orPair1=cell(1,NN);
% % maxMP = zeros(1,NN);
% % % pvals =zeros(1,NN);
% % % tic
% % for k=1:NN
% % %     k
% %     % win
% %     if ispc
% %     command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
% %        names2{k} ' --input_b_file_name=.\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
% %        ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
% %        ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
% %     else
% %        command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
% %        names2{k} ' --input_b_file_name=\' namesBar{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
% %        ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
% %        '