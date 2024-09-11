function [ZZ,data] = fit_all_constants(kN,kY, cN,cY, sigma, ligandLength,cI,bI,fastaFileF,parl,sets)

%% SQP?
[yoyoBindingProb] = get_yoyo_prob(kN,kY, cN,cY, sigma, ligandLength);

tol = 1e-11;
options_all = optimoptions(@fmincon,'Display', 'iter', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient',false, 'CheckGradient', false,  'OptimalityTolerance', tol, 'MaxFunctionEvaluations',10000 , 'StepTolerance', tol);

xopt = [yoyoBindingProb];

nvars = length(xopt);


import Core.rescale_barcode_data; % re-scale initial data
for i=1:length(bI)
    [barGenRe{i}] = rescale_barcode_data(bI{i},1,cI{i}.bestBarStretch);
end

data = cell(1,length(bI));
for i=1:length(bI)
    [~,mid,~] = fileparts(fastaFileF{i});
    data{i} = load(['seq_example',mid,'_',num2str(parl{i}{11}),'_',num2str(parl{i}{2}) ,'.mat']); % add nmpx att the end
end




 % yoyoBindingProb,compI,barGen,dataName,pars,sets
fun = @(xSv,ix1) -constfun_fit(xSv, cI{ix1},barGenRe{ix1},data{ix1},parl{ix1},sets);

funM = @(x) -sqrt(sum(arrayfun(@(i) fun(x,i),1:24).^2));

%%
% no bounds ...
lbounds = [zeros(1, nvars)];
ubounds = [0.3*ones(1, nvars)]; % no upper bounds on amplitudes


[ZZ] = fmincon(@(x) funM(x), xopt, [], [], [], [], lbounds, ubounds, [], options_all);

% figure,plot([ZZ-xopt]);hold on;plot(xopt)
%         [sortedSubseq, sortv, countATs, orderSeq,sortord] = sorted_NT(ligandLength);
% 
%         figure,plot(sortv*max(ZZ)/4);hold on
%         plot(ZZ(fliplr(sortord)))
% 
%         ZZ


end

