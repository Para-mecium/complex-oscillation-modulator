clear
clc

scriptDir = fileparts(mfilename('fullpath'));
resultsRootDir = fullfile(scriptDir, 'results');
scaleLevels = [1];

paramNames = {'R_C', 'C_1', 'C_2', 'C_3'};
xLabels = {'{\itR_C} (k\Omega)', '{\itC}_1 (\muF)', '{\itC}_2 (\muF)', '{\itC}_3 (\muF)'};

for scaleIdx = 1:numel(scaleLevels)
    scale = scaleLevels(scaleIdx);
    resultScaleDir = fullfile(resultsRootDir, scale_dir_name(scale));
    statsFile = fullfile(resultScaleDir, 'successful_params_stats.mat');

    data = load(statsFile);
    params = data.successPhysicalParams;
    corrMatrix = data.corrMatrixFiltered;
    statsSummary = data.statsSummary;

%% Histograms
histogramBinCount = 50;
histogramFaceColor = [43, 96, 166] / 255;
histogramEdgeColor = [6, 16, 111] / 255;
figHist = figure('Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:4
    ax = nexttile;
    values = params(:, i);
    histogram(ax, values, histogramBinCount, ...
        'Normalization', 'probability', ...
        'FaceColor', histogramFaceColor, ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', histogramEdgeColor, ...
        'EdgeAlpha', 1, ...
        'LineWidth', 0.5);
    hold(ax, 'on')
    grid(ax, 'on')
    box(ax, 'on')

    switch i
        case 1
            statsCurr = statsSummary.R_C;
        case 2
            statsCurr = statsSummary.C_1;
        case 3
            statsCurr = statsSummary.C_2;
        otherwise
            statsCurr = statsSummary.C_3;
    end

    % xline(ax, statsCurr.mean, '-', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
    % xline(ax, statsCurr.median, '--', 'LineWidth', 1.2, 'Color', [0.00 0.45 0.74]);

    ylim(ax, [0, 0.6])
    set(ax, 'YTick', 0:0.2:0.6)
    xlabel(ax, xLabels{i}, 'FontName', 'Arial')
    ylabel(ax, 'Probability', 'FontName', 'Arial')
    set(ax, 'FontName', 'Arial', 'FontSize', 10)
end

sgtitle(figHist, sprintf('Sensitivity to initial data, sampling region scale = %.3g', scale), ...
    'FontName', 'Arial', 'FontWeight', 'normal');

exportgraphics(figHist, fullfile(resultScaleDir, 'successful_params_histograms.png'), 'Resolution', 300);

%% Boxplot
figBox = figure('Color', 'w');
boxplot(params, 'Labels', paramNames, 'Symbol', 'k+');
grid on
box on
ylabel('Value', 'FontName', 'Arial')
set(gca, 'FontSize', 10)
title(sprintf('Sampling region scale = %.3g', scale), ...
    'FontName', 'Arial', 'FontWeight', 'normal');
exportgraphics(figBox, fullfile(resultScaleDir, 'successful_params_boxplot.png'), 'Resolution', 300);

%% Correlation heatmap
% figCorr = figure('Color', 'w');
% imagesc(corrMatrix);
% axis square
% colorbar
% colormap(parula)
% clim([-1 1])
% set(gca, 'XTick', 1:4, 'XTickLabel', paramNames, ...
%     'YTick', 1:4, 'YTickLabel', paramNames, 'FontSize', 10)
% 
% for i = 1:4
%     for j = 1:4
%         textColor = 'w';
%         if corrMatrix(i, j) > 0
%             textColor = 'k';
%         end
%         text(j, i, sprintf('%.2f', corrMatrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'Color', textColor, 'FontSize', 10);
%     end
% end
% 
% exportgraphics(figCorr, fullfile(resultScaleDir, 'successful_params_correlation_filtered.png'), 'Resolution', 300);

%% Scatter views in the physical parameter space
figScatter = figure('Color', 'w');

leftAx = axes(figScatter, 'Position', [0.06, 0.13, 0.5, 0.7]);
plot3(leftAx, params(:, 2), params(:, 3), params(:, 1), '.', ...
    'MarkerSize', 12, 'Color', [0.00 0.45 0.74]);
grid on
box on
xlabel('{\itC}_1 (\muF)', 'FontName', 'Arial')
ylabel('{\itC}_2 (\muF)', 'FontName', 'Arial')
zlabel('{\itR_C} (\Omega)', 'FontName', 'Arial')
set(leftAx, 'FontSize', 10)
view(-70, 40);

rightAx = axes(figScatter, 'Position', [0.65, 0.4, 0.3, 0.32]);
plot(rightAx, params(:, 4), params(:, 1), '.', ...
    'MarkerSize', 8, 'Color', [0.85 0.33 0.10]);
grid on
box on
xlabel('{\itC}_3 (\muF)', 'FontName', 'Arial')
ylabel('{\itR_C} (k\Omega)', 'FontName', 'Arial')
set(rightAx, 'FontSize', 8)

exportgraphics(figScatter, fullfile(resultScaleDir, 'successful_params_scatter.png'), 'Resolution', 300);

fprintf('scale=%.6g: exported parameter distribution figures to %s\n', scale, resultScaleDir);
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end
