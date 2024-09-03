function [ZZ] = run_single_parameter(gcSF,isC,sigma,kN,psf,cY,cN,kY,ligandLength,fastaFileF,expNr,...
    expPar,bI,cI,bound,sets)



nmbp = expPar{expNr}.nmbp;
nmpx = expPar{expNr}.nmpx;
[hcaSets,timestamp] = set_def();
hcaSets.pixelWidthNm = nmpx;
fastaFile = fastaFileF{expNr};
hcaSets.folder{1} = fastaFileF{expNr} ;

pxSize = nmpx/nmbp; % pixel size

parlist = @(psf) [gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY,cN,kY,ligandLength];
parlistcell = @(psf) num2cell(parlist(psf)) ;

%   [~,mid,~] = fileparts(hcaSets.folder{1} );
% delete(['seq_example',mid,'_',num2str(ligandLength) ,'.mat']);


% initial yoyo binding prob
[yoyoBindingProb] = get_yoyo_prob(kN,kY, cN,cY, sigma, ligandLength);
yoyofun = @(sigma)  get_yoyo_prob(kN,kY, cN,cY, sigma, ligandLength);



import Core.rescale_barcode_data; % re-scale initial data
[barGenRe] = rescale_barcode_data(bI{expNr},1,cI{expNr}.bestBarStretch);

fun = @(sx) -constfun_fit(yoyofun(sx), cI{expNr},barGenRe,hcaSets.folder{1},parlistcell(psf),sets);


tol = 1e-11;
options_all = optimoptions(@fmincon,'Display', 'iter', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient',false, 'CheckGradient', false,  'OptimalityTolerance', tol, 'MaxFunctionEvaluations',20000 , 'StepTolerance', tol);

xopt = [sigma];

lbounds = [bound(1)];
ubounds = [bound(2)]; % no upper bounds on amplitudes

[ZZ,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) fun(x), xopt, [], [], [], [], lbounds, ubounds, [], options_all);

end

