clear
clc

scriptDir = fileparts(mfilename('fullpath'));
statsFile = fullfile(scriptDir, 'successful_params_stats.mat');

data = load(statsFile);
params = data.successPhysicalParams;
corrMatrix = data.corrMatrix;
statsSummary = data.statsSummary;

paramNames = {'R_C', 'C_1', 'C_2', 'C_3'};
xLabels = {'{\itR_C} (k\Omega)', '{\itC}_1 (\muF)', '{\itC}_2 (\muF)', '{\itC}_3 (\muF)'};

%% Histograms
figHist = figure('Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:4
    nexttile
    histogram(params(:, i), 10, 'FaceColor', [0.72 0.80 0.90], 'EdgeColor', [0.25 0.25 0.25]);
    hold on
    grid on
    box on

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

    xline(statsCurr.mean, '-', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
    xline(statsCurr.median, '--', 'LineWidth', 1.2, 'Color', [0.00 0.45 0.74]);

    xlabel(xLabels{i}, 'FontName', 'Arial')
    ylabel('Count', 'FontName', 'Arial')
    title(paramNames{i}, 'FontName', 'Arial')
    set(gca, 'FontSize', 10)
end

exportgraphics(figHist, fullfile(scriptDir, 'successful_params_histograms.png'), 'Resolution', 300);

%% Boxplot
figBox = figure('Color', 'w');
boxplot(params, 'Labels', paramNames, 'Symbol', 'k+');
grid on
box on
ylabel('Value', 'FontName', 'Arial')
set(gca, 'FontSize', 10)
exportgraphics(figBox, fullfile(scriptDir, 'successful_params_boxplot.png'), 'Resolution', 300);

%% Correlation heatmap
figCorr = figure('Color', 'w');
imagesc(corrMatrix);
axis square
colorbar
colormap(parula)
clim([-1 1])
set(gca, 'XTick', 1:4, 'XTickLabel', paramNames, ...
    'YTick', 1:4, 'YTickLabel', paramNames, 'FontSize', 10)
title('Correlation Matrix', 'FontName', 'Arial')

for i = 1:4
    for j = 1:4
        text(j, i, sprintf('%.2f', corrMatrix(i, j)), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 10);
    end
end

exportgraphics(figCorr, fullfile(scriptDir, 'successful_params_correlation.png'), 'Resolution', 300);

%% Scatter views in the physical parameter space
figScatter = figure('Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile
plot3(params(:, 2), params(:, 3), params(:, 1), '.', 'MarkerSize', 12, 'Color', [0.00 0.45 0.74]);
grid on
box on
xlabel('{\itC}_1 (\muF)', 'FontName', 'Arial')
ylabel('{\itC}_2 (\muF)', 'FontName', 'Arial')
zlabel('{\itR_C} (k\Omega)', 'FontName', 'Arial')
set(gca, 'FontSize', 10)
view(-70, 40);

nexttile
plot(params(:, 4), params(:, 1), '.', 'MarkerSize', 12, 'Color', [0.85 0.33 0.10]);
grid on
box on
xlabel('{\itC}_3 (\muF)', 'FontName', 'Arial')
ylabel('{\itR_C} (k\Omega)', 'FontName', 'Arial')
set(gca, 'FontSize', 10)

exportgraphics(figScatter, fullfile(scriptDir, 'successful_params_scatter.png'), 'Resolution', 300);
