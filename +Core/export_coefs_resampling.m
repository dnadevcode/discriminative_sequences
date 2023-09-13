function [T] = export_coefs_resampling(scores, barGen, mpval,matDirpath,timestamp)

    %
    %   Args:
    %       scores, barGen, mpval,matDirpath

    % Returns:
    %   T - table

    % export_coefs_resampling
    %     w - window width
    % bootstrapped PCC - average over alignment scores score_b for N=1000 bootstrap/concatenated fragments of length 20 px,
    % s_b - sigma over alignment scores score_b
    % N_T - number of theories within 0.05 of max PCC
    % N_{T,Unique} - number of unique theories withiin 0.05 of max PCC
    % ms - match score (-norminv(pval)), where pval=1-evcdf(maxpcc,nu*lenOverlap, 2*max(shortlength,longlength)). nu is somewhere between 0.04-0.06.
    % ms_b - average over bootstrapped match scores
    % s_{ms_b} - sigma over bootstrapped match scores

    rows = arrayfun(@(x) ['w=',num2str(x)], mpval,'un',false);
    bar = arrayfun(@(x) matlab.lang.makeValidName(strrep(barGen{1}.name,'.tif','')),1:length(rows),'un',false);


    T1 = table(bar',mpval','VariableNames',{'Intensity profile','w'});
    T2 = array2table(scores,'VariableNames',{'bootstrapped PCC','s_b','N_T', 'N_{T,unique}','ms','ms_b','s_{ms_b}'});

    T = [T1 T2];

    writemode = 'append';
    disp('Saving ccvals table');
    import CBT.Hca.Export.export_cc;
    export_cc(T, matDirpath,timestamp,writemode);      


end

