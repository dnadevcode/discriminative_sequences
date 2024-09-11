function [ZZ] = fit_two_constants(kN,kY, cN,cY, sigma, ligandLength,cI,bI,fastaFileF,parl,sets)
% fit sigma and psf:

%% SQP?
[yoyoBindingProb] = get_yoyo_prob(kN,kY, cN,cY, sigma, ligandLength);
yoyofun = @(sigma)  get_yoyo_prob(kN, kY, cN, cY, sigma, ligandLength);

tol = 1e-11;
options_all = optimoptions(@fmincon,'Display', 'iter', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient',false, 'CheckGradient', false,  'OptimalityTolerance', tol, 'MaxFunctionEvaluations',1000 , 'StepTolerance', tol);

xopt = [yoyoBindingProb];

nvars = length(xopt);

parnr = 7;


import Core.rescale_barcode_data; % re-scale initial data
for i=1:length(bI)
    [barGenRe{i}] = rescale_barcode_data(bI{i},1,cI{i}.bestBarStretch);
end

% parlist =  [gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY,cN,kY,ligandLength];
parlistcell = @(x, parnr, ix) [parl{ix}(1:parnr-1) x parl{ix}(parnr+1:end)] ;

data = cell(1,length(bI));
for i=1:length(bI)
    [~,mid,~] = fileparts(fastaFileF{i});
    data{i} = load(['seq_example',mid,'_',num2str(parl{i}{11}),'_',num2str(parl{i}{2}) ,'.mat']); % add nmpx att the end
end

 % yoyoBindingProb,compI,barGen,dataName,pars,sets
% fun = @(xSv,ix1) -constfun_fit(xSv, cI{ix1},bI{ix1},fastaFileF{ix1},parl{ix1},sets);
fun = @(sigma, x, ix1) -constfun_fit(yoyofun(sigma), cI{ix1}, barGenRe{ix1}, data{ix1}, parlistcell(x, parnr,ix1), sets);

funM = @(sigma,x) -sqrt(sum(arrayfun(@(i) fun(sigma,x,i),1:24).^2));


%%
% no bounds ...
lbounds = [0 100];
ubounds = [5 600]; % no upper bounds on amplitudes

xopt = [0.4, 370];

[ZZ] = fmincon(@(x) funM(x(1), x(2)), xopt, [], [], [], [], lbounds, ubounds, [], options_all);

% figure,plot([ZZ-xopt]);hold on;plot(xopt)
%         [sortedSubseq, sortv, countATs, orderSeq,sortord] = sorted_NT(ligandLength);
% 
%         figure,plot(sortv*max(ZZ)/4);hold on
%         plot(ZZ(fliplr(sortord)))
% 
%         ZZ


end

