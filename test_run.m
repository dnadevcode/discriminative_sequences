% example run:
ix = 1;
iy = 1;
dirName = '/export/scratch/albertas/download_dump/S. pyogenes all data/onlyKymo/';
% dirName = '/export/scratch/albertas/download_dump/pyo_test/';
depth=1
windowWidths=[0 250:50:1000];
sF = 0.9:0.025:1.1;
thryFiles = dir('/export/scratch/albertas/download_dump/single/theoryOutput/*.mat');


local_pipeline_mp(ix, iy, dirName, depth, windowWidths, sF, thryFiles)

 c = parcluster;
% c.AdditionalProperties.AccountName = 'snic2021-5-132';
c.AdditionalProperties.WallTime = '6:00:00';
% j =batch(c,@local_pipeline_mp,1,{ix, iy, dirName, depth, 0, sF, thryFiles},'Pool',30);

%%
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/hca'));
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/discriminative_sequences'));



dirName = '/proj/snic2022-5-384/users/x_albdv/data/discriminative/pyo/';

dr = dir(dirName);
folderNames = {dr([dr.isdir]).name};
dr = dr(~ismember(folderNames ,{'.','..'}));

numDirs = length(dr);
subdirs = zeros(1,numDirs);
barN = cell(1,numDirs);
for i=1:numDirs
    dr2 = dir(fullfile(dr(i).folder,dr(i).name));
    dr2 = dr2([dr2.isdir]);
    folderNames = {dr2([dr2.isdir]).name};
    dr2 = dr2(~ismember(folderNames ,{'.','..'}));
    subdirs(i) = numel(dr2);
    for j=1:numel(dr2)
        barN{i}(j) = numel(dir(fullfile(dr2(j).folder,dr2(j).name,'kymos','*.tif')));
    end
end


twoList = zeros(sum(subdirs),2);
id=1;
for i=1:length(subdirs)
    for j=1:subdirs(i)
        twoList(id,:) = [i j];
        id = id+1;
    end
end


% total sum(subdirs) folders to run. First run everything at 0 to check
windowWidths=0 ; %[0 250:50:1000];
sF = 0.9:0.025:1.1;
thryFiles = dir('/proj/snic2022-5-384/users/x_albdv/data/CHR/New Ref-Theory files April 2024/*.mat');
depth=1;

v = cell(numDirs,subdirs(i));
for i=4:numDirs
    for j=1:subdirs(i)
        v{i,j} =batch(c,@local_pipeline_mp,1,{i, j, dirName, depth, 0, sF, thryFiles},'Pool',30);
    end
end

v = cell(numDirs,subdirs(i));
w=400;
for i=1:numDirs
    for j=1:subdirs(i)
        v{i,j} =batch(c,@local_pipeline_mp,1,{i, j, dirName, depth, w, sF, thryFiles},'Pool',60);
    end
end

%%
% example run:
ix=1
iy=1

 c = parcluster;
% c.AdditionalProperties.AccountName = 'snic2021-5-132';
% c.AdditionalProperties.WallTime = '1:00:00';


local_pipeline_mp(ix, iy, dirName, depth, windowWidths, sF, thryFiles)

% re-calculate the ones that failed before
w = 350;
batch(c,@local_pipeline_mp,1,{ix,iy, dirName, depth, w, sF, thryFiles},'Pool',29);

% re-calculate the ones that failed before
w = 300;
idd = 22;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',29);

%
w = 450;
v = cell(numDirs,subdirs(i));
for i=1:numDirs
    for j=1:subdirs(i)
        v{i,j} =batch(c,@local_pipeline_mp,1,{i, j, dirName, depth, w, sF, thryFiles},'Pool',29);
    end
end


%%
numW = 60;
w = [750 800];
idd = 1;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);

w = [800:50:1000];
idd = 2;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);


w = [ 750 800];
idd = 3;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [  750 800];
idd = 4;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);


w = [ 800];
idd = 5;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [ 750 800];
idd = 6;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [ 450 550 750 800];
idd = 7;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [250 450 500 750 800];
idd = 8;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);


w = [ 450  750 800];
idd = 9;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);





w = [ 550 750 800];
idd = 10;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);





w = [  500 750 800];
idd = 11;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);






w = [    750 800];
idd = 12;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);






w = [  450 500 650:50:1000];
idd = 13;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);




w = [  550  750 800];
idd = 14;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);

%

w = [  450 550  750 800];
idd = 15;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);


w = [  450  550 750 800 900:50:100];
idd = 16;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [   750:50:1000];
idd = 17;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);


%

w = [   450 500 550 750 800];
idd = 18;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);





w = [ 250  750 800];
idd = 19;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);





%


w = [  500 750 800];
idd = 20;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);


%%
w = [ 250 500 750 800];
idd = 21;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [ 250 300 450 500 750 800];
idd = 22;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [ 450:50:1000];
idd = 23;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);



w = [ 500:50:1000];
idd = 24;
batch(c,@local_pipeline_mp,1,{twoList(idd,1),twoList(idd,2), dirName, depth, w, sF, thryFiles},'Pool',numW);




