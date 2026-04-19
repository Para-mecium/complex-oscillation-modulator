clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'data', 'fig3f', 'fig3f_data.mat');
figureFile = fullfile(scriptDir, 'fig3f.png');

%% Load data
loaded = load(dataFile);
startPoint = loaded.startPoint;
orthPeriod = loaded.orthPeriod;
orthAmplitude = loaded.orthAmplitude;
directMid = loaded.directMid;
directPath = loaded.directPath;

%% Build path coordinates
directMidParams = [startPoint.Parameters; vertcat(directMid.logs.params)];
if size(directMidParams, 1) >= 2 && norm(directMidParams(2, :) - startPoint.Parameters, inf) < 1e-10
    directMidParams(2, :) = [];
end

directEndParams = [directMid.Parameters; vertcat(directPath.logs.params)];
if size(directEndParams, 1) >= 2 && norm(directEndParams(2, :) - directMid.Parameters, inf) < 1e-10
    directEndParams(2, :) = [];
end

directPathParams = [directMidParams; directEndParams];
if size(directPathParams, 1) >= 2 && norm(directPathParams(size(directMidParams, 1) + 1, :) - directMidParams(end, :), inf) < 1e-10
    directPathParams(size(directMidParams, 1) + 1, :) = [];
end

orthPeriodParams = [startPoint.Parameters; vertcat(orthPeriod.logs.params)];
if size(orthPeriodParams, 1) >= 2 && norm(orthPeriodParams(2, :) - startPoint.Parameters, inf) < 1e-10
    orthPeriodParams(2, :) = [];
end

orthAmplitudeParams = [orthPeriod.Parameters; vertcat(orthAmplitude.logs.params)];
if size(orthAmplitudeParams, 1) >= 2 && norm(orthAmplitudeParams(2, :) - orthPeriod.Parameters, inf) < 1e-10
    orthAmplitudeParams(2, :) = [];
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 80, 1380, 760]);
tiled = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

axPath = nexttile(tiled, [3, 1]);
hold(axPath, 'on');
plot(axPath, directPathParams(:, 1), directPathParams(:, 2), ...
    'Color', [0.42, 0.64, 0.1], 'LineWidth', 3, 'DisplayName', 'Direct modulation path');
plot(axPath, orthPeriodParams(:, 1), orthPeriodParams(:, 2), ...
    'Color', [0.98, 0.72, 0.25], 'LineWidth', 3, 'DisplayName', 'Orthogonal period path');
plot(axPath, orthAmplitudeParams(:, 1), orthAmplitudeParams(:, 2), ...
    'Color', [0.42, 0.7, 0.92], 'LineWidth', 3, 'DisplayName', 'Orthogonal amplitude path');
scatter(axPath, startPoint.Parameters(1), startPoint.Parameters(2), 220, [0.12, 0.5, 0.85], 'filled');
scatter(axPath, directMid.Parameters(1), directMid.Parameters(2), 220, [0.85, 0.85, 0.85], 'filled', 'd');
scatter(axPath, orthPeriod.Parameters(1), orthPeriod.Parameters(2), 220, [0.85, 0.85, 0.85], 'filled', 's');
scatter(axPath, directPath.Parameters(1), directPath.Parameters(2), 320, [0.98, 0.78, 0.18], 'filled', 'p');
grid(axPath, 'on');
xlabel(axPath, 'E_T (a.u.)');
ylabel(axPath, 'Temperature (K)');
legend(axPath, 'Location', 'northwest');
title(axPath, 'Direct vs orthogonal implementations');

plotData = cell(1, 6);
plotData{1} = startPoint;
plotData{2} = startPoint;
plotData{3} = directMid;
plotData{4} = orthPeriod;
plotData{5} = directPath;
plotData{6} = orthAmplitude;

plotColors = cell(1, 6);
plotColors{1} = [0.12, 0.5, 0.85];
plotColors{2} = [0.12, 0.5, 0.85];
plotColors{3} = [0.86, 0.86, 0.86];
plotColors{4} = [0.86, 0.86, 0.86];
plotColors{5} = [0.98, 0.78, 0.18];
plotColors{6} = [0.98, 0.78, 0.18];

plotTitles = cell(1, 6);
plotTitles{1} = 'Direct start';
plotTitles{2} = 'Orthogonal start';
plotTitles{3} = 'Direct middle';
plotTitles{4} = 'Orthogonal middle';
plotTitles{5} = 'Direct end';
plotTitles{6} = 'Orthogonal end';

for i = 1:6
    ax = nexttile(tiled);
    hold(ax, 'on');

    t = plotData{i}.TS{1}(:);
    y = plotData{i}.TS{2};
    tRel = t - t(1);
    tStep = tRel(end);
    tLong = tRel;
    yLong = y;
    for k = 1:5
        tLong = [tLong; tRel(2:end) + k * tStep]; %#ok<AGROW>
        yLong = [yLong; y(2:end, :)]; %#ok<AGROW>
    end
    idxMax = find(islocalmax(yLong(:, 2)));
    i0 = idxMax(end - 1);
    i1 = idxMax(end);
    tPlot = tLong(i0:i1) - tLong(i0);
    yPlot = yLong(i0:i1, 2);

    plot(ax, tPlot, yPlot, 'LineWidth', 3, 'Color', plotColors{i});
    grid(ax, 'on');
    xlabel(ax, 'Time (min)');
    ylabel(ax, 'Y (a.u.)');
    title(ax, plotTitles{i});
end

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);
