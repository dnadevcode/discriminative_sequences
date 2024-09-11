newDir = '/export/scratch/albertas/output_dump/outputPCCGroupedEHEC/';
mkdir(newDir);
merge_alignment_outputs(saveDir,newDir, twoList)


addpath(genpath('/proj/dnadevdata/reps/discriminative_sequences'))
cdiff = 0.05; % c diff
sflevel = 1; % 1 :-10.. 2: -7.5.. 3: -5 .. ...
[distinctoutput, mostCommonSequence,mostCommonRep,namesDir,thryDiscName] = get_best_theory(newDir, [], idSpecies(thryId),barN, twoList,foldname,thryNames(thryId),'hca', sflevel,cdiff);

txt = '/home/avesta/albertas/reps/discriminative_sequences/constants/ehec_t.txt';

idxs = cellfun(@(y) str2num(y{1}),cellfun(@(x) strsplit(x,'_'),namesDir,'UniformOutput',false));

[a,sortId] = sort(idxs);

fd = fopen(txt,'w');
cellfun(@(x) fprintf(fd,'%s\n',lower(x)),thryDiscName(sortId));
fclose(fd);

%% validation
% dirNameV = '/proj/dnadevdata/users/x_albdv/data/EHEC/EHEC data for local alignment/';
dirNameV = '/export/scratch/albertas/data_temp/Alignment/data/ehec/EHEC data for local alignment/';
refsNamesV = 'ehec_t.txt'
[kymoStructsV,barNV,twoListV,bGV,expParV,fastaFileFV] = load_kymo_data_from_fold(dirNameV, refsNamesV,allTheoryFold,0.8:0.025:1.2,0);

%%
twoParameterFit = ZZ;

[cIV,bIV,compI2V, parlV,allCoefsFitV,m2V,mAV] = get_pseudotheory_positions(bGV, twoListV, expParV, fastaFileFV, sets, twoParameterFit(2), twoParameterFit(1));

mtwoV = cellfun(@(x) mean(x.maxcoef),cIV,'UniformOutput',true);

plot_at_preference_model(bGV,mAV,mtwoV,'EHEC')