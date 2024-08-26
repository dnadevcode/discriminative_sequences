function [m] = constfun_fit(yoyoBindingProb,compI,barGen,dataName,pars,sets)

    import Helper.get_theory_twostate_fit;

%     k = 1;


import CBT.Hca.Core.Comparison.hca_compare_distance;
import Core.Discriminative.generate_sf_struct;

    
    [~,thr] = get_theory_twostate_fit(dataName,pars{:},yoyoBindingProb);

    if ~isempty(compI)
%         tic
        compInew = [];
%         tic
        for i=1:length(barGen)
            theory = thr{1};
            theory.rawBarcode = [ theory.rawBarcode  theory.rawBarcode(1:length(barGen{i}.rescaled{1}.rawBarcode))];
            theory.rawBarcode = theory.rawBarcode(compI.pos(i):compI.pos(i)+length(barGen{i}.rescaled{1}.rawBarcode)-1);
            theory.length = length( theory.rawBarcode);
% 
            a = theory.rawBarcode(barGen{i}.rescaled{1}.rawBitmask);
            b = barGen{i}.rescaled{1}.rawBarcode(barGen{i}.rescaled{1}.rawBitmask);
            if compI.or(i) ==2
                b= fliplr(b);
            end
            compInew{i}.maxcoef = zscore(a,1)*zscore(b',1)/length(a);


            % todo - length rescaleand directly compare

%             [rezMax] = hca_compare_distance(barGen(i),theory, sets );
% 
%             [compInew{i}] = generate_sf_struct(rezMax,sets);
%     

%             [~,rez{i},~] = compare_to_t(barGen(i),theory,cp{i}.bestBarStretch,sets);
        end
         mvals = cellfun(@(y) y.maxcoef, compInew);

%     toc
%             m
% %         tomc
        
    else
%         tic
% tic
        [rezMax] = hca_compare_distance(barGen,thr{1}, sets );

        [compInew] = generate_sf_struct(rezMax,sets);
        mvals = compInew.maxcoef;
% toc
%         [~,rezI{k},~] = compare_to_t(barGen,theoryStr{k},sF,sets);
%         toc
%         m = cellfun(@(y) mean(cellfun(@(x) x.maxcoef(1),y)), rezI{k});
%         m
    end
    m = mean(mvals);
% 
%     figure,plot( zscore(theory.rawBarcode))
%     hold on
%     plot(zscore(barGen{1}.rescaled{1}.rawBarcode))
%     a = theory.rawBarcode(barGen{1}.rescaled{1}.rawBitmask);
%     b = barGen{1}.rescaled{1}.rawBarcode(barGen{1}.rescaled{1}.rawBitmask);
%     zscore(a,1)*zscore(b',1)/length(a)

end

