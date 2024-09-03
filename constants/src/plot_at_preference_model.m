function plot_at_preference_model(bG,mA,mtwo)
    % Comparison plot
    % lengths = cellfun(@(x) length(x),bI);
    figure,bar([mA(1:length(bG));mtwo(1:length(bG))]')
    hold on
    title('Models comparison')
    legend({'hca','2state'},'location','southoutside')
    xlabel('S.Pyo')
    ylabel('Mean PCC')
    ylim([0.5 0.8])
end
