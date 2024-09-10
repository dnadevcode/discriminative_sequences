function plot_at_preference_model(bG,mA,mtwo,name)

    if nargin < 4
        name = 'S.Pyogenes';
    end
    % Comparison plot
    % lengths = cellfun(@(x) length(x),bI);
    figure,bar([mA(1:length(bG));mtwo(1:length(bG))]')
    hold on
    title('Models comparison')
    legend({'hca','2state'},'location','southoutside')
    xlabel(name)
    ylabel('Mean PCC')
    ylim([0.5 0.8])
end
