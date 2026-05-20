clear
clc

scriptDir = fileparts(mfilename('fullpath'));
resultsRootDir = fullfile(scriptDir, 'results');
scaleLevels = [0.5 0.75 1 1.25 1.5];

scaleValues = scaleLevels(:);
numTotal = zeros(numel(scaleValues), 1);
numSuccess = zeros(numel(scaleValues), 1);
successRates = zeros(numel(scaleValues), 1);

for scaleIdx = 1:numel(scaleValues)
    scale = scaleValues(scaleIdx);
    summaryFile = fullfile(resultsRootDir, scale_dir_name(scale), ...
        'params_inf_sensitivity_summary.mat');
    summaryData = load(summaryFile);

    numTotal(scaleIdx) = summaryData.numTotal;
    numSuccess(scaleIdx) = summaryData.numSuccess;
    successRates(scaleIdx) = summaryData.successRate;

    fprintf('scale=%.6g: continuation success=%d/%d (%.3f)\n', ...
        scale, numSuccess(scaleIdx), numTotal(scaleIdx), successRates(scaleIdx));
end

fig = figure();
ax = axes(fig);
hold(ax, 'on')
grid(ax, 'on')
box(ax, 'on')

lineColor = [0.00 0.45 0.74];
plot(ax, scaleValues, successRates, 'o-', ...
    'Color', lineColor, ...
    'MarkerFaceColor', lineColor, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.4, ...
    'MarkerSize', 6);

for scaleIdx = 1:numel(scaleValues)
    text(ax, scaleValues(scaleIdx), min(successRates(scaleIdx) + 0.04, 1.02), ...
        sprintf('%d/%d', numSuccess(scaleIdx), numTotal(scaleIdx)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', ...
        'FontSize', 9, ...
        'Color', [0.20 0.20 0.20]);
end

xlim(ax, [min(scaleValues) - 0.05, max(scaleValues) + 0.05]);
ylim(ax, [0.4, 1.08]);
set(ax, ...
    'XTick', scaleValues, ...
    'XTickLabel', compose('%.3g', scaleValues), ...
    'YTick', 0:0.2:1, ...
    'YTickLabel', compose('%d%%', 0:20:100), ...
    'FontName', 'Arial', ...
    'FontSize', 11);
xlabel(ax, 'Sampling region scale', 'FontName', 'Arial');
ylabel(ax, 'Modulation success rate', 'FontName', 'Arial');
title(ax, 'Sensitivity to initial data', ...
    'FontName', 'Arial', 'FontWeight', 'normal');

outputPng = fullfile(resultsRootDir, 'Sensitivity_to_initial_data_success_rate.png');
outputPdf = fullfile(resultsRootDir, 'Sensitivity_to_initial_data_success_rate.pdf');
exportgraphics(fig, outputPng, 'Resolution', 300);
exportgraphics(fig, outputPdf, 'ContentType', 'vector');

successRateTable = table(scaleValues, numSuccess, numTotal, successRates, ...
    'VariableNames', {'scale', 'numSuccess', 'numTotal', 'successRate'});
save(fullfile(resultsRootDir, 'Sensitivity_to_initial_data_success_rate.mat'), ...
    'successRateTable');

fprintf('Saved success-rate curve: %s\n', outputPng);
fprintf('Saved vector figure: %s\n', outputPdf);

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end
