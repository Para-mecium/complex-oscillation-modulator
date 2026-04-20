clear
clc

scriptDir = fileparts(mfilename('fullpath'));
statsFile = fullfile(scriptDir, 'forward_perturbation_stats.mat');

data = load(statsFile);
scenarioLabels = data.scenarioLabels;
periodRelChange = data.periodRelChange;
varAmpRelChange = data.varAmpRelChange;

%% Period change
figPeriod = figure('Color', 'w');
bar(100 * periodRelChange, 'FaceColor', [0.30 0.55 0.80], 'EdgeColor', [0.20 0.20 0.20]);
grid on
box on
yline(0, 'k-');
xlim([0.5, numel(scenarioLabels) + 0.5]);
set(gca, 'XTick', 1:numel(scenarioLabels), 'XTickLabel', scenarioLabels, ...
    'XTickLabelRotation', 45, 'FontSize', 10)
ylabel('Period Change (%)', 'FontName', 'Arial')
title('Relative Period Change Under Parameter Perturbations', 'FontName', 'Arial')
exportgraphics(figPeriod, fullfile(scriptDir, 'forward_perturb_period_change.png'), 'Resolution', 300);

%% Amplitude changes
figAmp = figure('Color', 'w');
bar(100 * varAmpRelChange, 'grouped');
grid on
box on
yline(0, 'k-');
xlim([0.5, numel(scenarioLabels) + 0.5]);
set(gca, 'XTick', 1:numel(scenarioLabels), 'XTickLabel', scenarioLabels, ...
    'XTickLabelRotation', 45, 'FontSize', 10)
ylabel('Amplitude Change (%)', 'FontName', 'Arial')
title('Relative Amplitude Change Under Parameter Perturbations', 'FontName', 'Arial')
legend({'V_1', 'V_2', 'V_3'}, 'Location', 'best')
exportgraphics(figAmp, fullfile(scriptDir, 'forward_perturb_varAmp_change.png'), 'Resolution', 300);

%% Absolute-change heatmap
figHeat = figure('Color', 'w');
h = heatmap(scenarioLabels, {'Period', 'Amp_1', 'Amp_2', 'Amp_3'}, ...
    100 * [abs(periodRelChange(:)).'; abs(varAmpRelChange).'], ...
    'Colormap', parula);
h.Title = 'Absolute Relative Change (%)';
exportgraphics(figHeat, fullfile(scriptDir, 'forward_perturb_heatmap.png'), 'Resolution', 300);
