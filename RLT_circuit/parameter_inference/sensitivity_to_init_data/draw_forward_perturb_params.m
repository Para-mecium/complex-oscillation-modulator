clear
clc

scriptDir = fileparts(mfilename('fullpath'));
statsFile = fullfile(scriptDir, 'forward_perturbation_stats.mat');

data = load(statsFile);
scenarioLabels = data.scenarioLabels;
periodRelChange = data.periodRelChange;
varAmpRelChange = data.varAmpRelChange;

%% Period and amplitude changes
figRelChange = figure('Color', 'w');
relChangePercent = 100 * [periodRelChange(:), varAmpRelChange];
bar(relChangePercent, 'grouped');
grid on
box on
yline(0, 'k-');
xlim([0.5, numel(scenarioLabels) + 0.5]);
set(gca, 'XTick', 1:numel(scenarioLabels), 'XTickLabel', scenarioLabels, ...
    'XTickLabelRotation', 45, 'FontSize', 10)
ylabel('Relative Change (%)', 'FontName', 'Arial')
title('Relative Period and Amplitude Change Under Parameter Perturbations', 'FontName', 'Arial')
legend({'Period', 'V_1 amplitude', 'V_2 amplitude', 'V_3 amplitude'}, 'Location', 'best')
exportgraphics(figRelChange, fullfile(scriptDir, 'forward_perturb_period_varAmp_change.png'), 'Resolution', 300);

%% Absolute-change heatmap
figHeat = figure('Color', 'w');
h = heatmap(scenarioLabels, {'Period', 'Amp_1', 'Amp_2', 'Amp_3'}, ...
    100 * [abs(periodRelChange(:)).'; abs(varAmpRelChange).'], ...
    'Colormap', parula);
h.Title = 'Magnitude of Relative Change (%)';
exportgraphics(figHeat, fullfile(scriptDir, 'forward_perturb_heatmap.png'), 'Resolution', 300);
