function [yoyoBindingProb] = get_yoyo_prob(kN,kY, cN,cY, sigma, ligandLength)


constFun = arrayfun(@(x) kN*exp(-x/sigma),(0:ligandLength));
probYoyo = cY*kY./(1+cY*kY+cN.*constFun);
[sortedSubseq, sortv, countATs, orderSeq] = sorted_NT(ligandLength);
yoyoBindingProb = ones(1,4^ligandLength);
yoyoBindingProb(orderSeq) =  probYoyo(ligandLength+1-countATs);
netropsinBindingConst = constFun(ligandLength+1-countATs);

   
end

