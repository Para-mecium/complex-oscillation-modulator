clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
figureFile = fullfile(scriptDir, 'fig3d.png');
bifurcationFile = fullfile(scriptDir, 'data', 'fig3', 'fig3_bifurcation_line.mat');
curveDataDir = fullfile(scriptDir, 'data', 'fig3d', 'curves');
markerDataDir = fullfile(scriptDir, 'data', 'fig3d', 'markers');

%% Files
targetAmplitudes = 1.2:0.3:3.6;
markerAmplitude = 2.1;
markerPeriods = [60, 80, 100];

curveFiles = cell(1, numel(targetAmplitudes));
for i = 1:numel(targetAmplitudes)
    amplitudeTag = strrep(sprintf('%.1f', targetAmplitudes(i)), '.', 'p');
    curveFiles{i} = fullfile(curveDataDir, sprintf('fig3d_iso_amplitude_A%s.mat', amplitudeTag));
end

markerFiles = cell(1, numel(markerPeriods));
markerAmplitudeTag = strrep(sprintf('%.1f', markerAmplitude), '.', 'p');
for i = 1:numel(markerPeriods)
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('fig3d_marker_A%s_T%03d.mat', markerAmplitudeTag, round(markerPeriods(i))));
end

%% Plot settings
lineColor = [0.25, 0.25, 0.25];
lineStyle = '--';
lineWidth = 1.4;

nColor = 256;
tColor = linspace(0, 1, nColor).';
coolColormap = (1 - tColor) .* [0.05, 0.25, 0.55] + tColor .* [0.88, 0.93, 0.97];

%% Load bifurcation lines
loaded = load(bifurcationFile);
visibleHopfCurve = loaded.visibleHopfCurve;
visibleLpCurve = loaded.visibleLpCurve;
xLimits = loaded.xRange;
yLimits = loaded.yRange;

hopfI1 = visibleHopfCurve.I1(:);
hopfET = visibleHopfCurve.ET(:);
hopfDx = abs(diff(hopfI1));
hopfDy = abs(diff(hopfET));
hopfDxRef = max(median(hopfDx(hopfDx > 0)), eps);
hopfDyRef = max(median(hopfDy(hopfDy > 0)), eps);
hopfBreakMask = hopfDx > 3 * hopfDxRef | hopfDy > max(0.12, 6 * hopfDyRef);
hopfI1Plot = hopfI1;
hopfETPlot = hopfET;
hopfI1Plot(find(hopfBreakMask) + 1) = NaN;
hopfETPlot(find(hopfBreakMask) + 1) = NaN;

lpI1 = visibleLpCurve.I1(:);
lpET = visibleLpCurve.ET(:);
lpDx = abs(diff(lpI1));
lpDy = abs(diff(lpET));
lpDxRef = max(median(lpDx(lpDx > 0)), eps);
lpDyRef = max(median(lpDy(lpDy > 0)), eps);
lpBreakMask = lpDx > 3 * lpDxRef | lpDy > max(0.12, 6 * lpDyRef);
lpI1Plot = lpI1;
lpETPlot = lpET;
lpI1Plot(find(lpBreakMask) + 1) = NaN;
lpETPlot(find(lpBreakMask) + 1) = NaN;

%% Load iso-amplitude curves
curveI1 = cell(1, numel(curveFiles));
curveET = cell(1, numel(curveFiles));
curvePeriod = cell(1, numel(curveFiles));
allCurvePeriods = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});

    seedParams = loaded.seed.params(:).';
    seedPeriod = loaded.seed.derivedView.period;

    leftParams = vertcat(loaded.leftBranch.params);
    leftPeriod = arrayfun(@(entry) entry.derived.period, loaded.leftBranch).';
    rightParams = vertcat(loaded.rightBranch.params);
    rightPeriod = arrayfun(@(entry) entry.derived.period, loaded.rightBranch).';

    curveI1{i} = [flipud(leftParams(:, 1)); seedParams(1); rightParams(:, 1)];
    curveET{i} = [flipud(leftParams(:, 2)); seedParams(2); rightParams(:, 2)];
    curvePeriod{i} = [flipud(leftPeriod); seedPeriod; rightPeriod];

    allCurvePeriods = [allCurvePeriods; curvePeriod{i}]; %#ok<AGROW>
end

periodLimits = [min(allCurvePeriods), max(allCurvePeriods)];
if periodLimits(1) >= periodLimits(2)
    delta = max(abs(periodLimits(1)) * 0.05, 1e-6);
    periodLimits = [periodLimits(1) - delta, periodLimits(2) + delta];
end

%% Load marker data
markerParams = cell(1, numel(markerFiles));
markerT = cell(1, numel(markerFiles));
markerY = cell(1, numel(markerFiles));
markerPeriodValues = zeros(1, numel(markerFiles));

for i = 1:numel(markerFiles)
    loaded = load(markerFiles{i});
    markerParams{i} = loaded.Parameters(:).';
    markerT{i} = loaded.TS{1};
    markerY{i} = loaded.TS{2};
    markerPeriodValues(i) = loaded.period;
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 1280, 520]);
tiled = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled, 1);
hold(ax1, 'on');
colormap(ax1, coolColormap);

for i = 1:numel(curveFiles)
    surface(ax1, [curveI1{i}.'; curveI1{i}.'], [curveET{i}.'; curveET{i}.'], ...
        zeros(2, numel(curveI1{i})), [curvePeriod{i}.'; curvePeriod{i}.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
end

for i = 1:numel(markerFiles)
    scatter(ax1, markerParams{i}(1), markerParams{i}(2), 170, markerPeriodValues(i), 'filled', ...
        'MarkerEdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.2);
end

plot(ax1, hopfI1Plot, hopfETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
plot(ax1, lpI1Plot, lpETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);

grid(ax1, 'on');
clim(ax1, periodLimits);
cb = colorbar(ax1);
cb.Label.String = 'Period (min)';
xlabel(ax1, 'I_1 (a.u.)');
ylabel(ax1, 'E_T (a.u.)');
title(ax1, sprintf('Iso-amplitude curves (A = %.1f highlighted)', markerAmplitude));
xlim(ax1, xLimits);
ylim(ax1, yLimits);

rightTiled = tiledlayout(tiled, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
rightTiled.Layout.Tile = 2;

plotDuration = 300;
xTicks = 0:60:plotDuration;
allMarkerY = vertcat(markerY{:});
concentrationLimits = [min(allMarkerY(:, 2)), max(allMarkerY(:, 2))];
concentrationPadding = 0.05 * diff(concentrationLimits);
concentrationLimits = concentrationLimits + [-concentrationPadding, concentrationPadding];
plotOrder = numel(markerFiles):-1:1;

for row = 1:numel(plotOrder)
    i = plotOrder(row);
    alpha = (markerPeriodValues(i) - periodLimits(1)) / (periodLimits(2) - periodLimits(1));
    alpha = min(max(alpha, 0), 1);
    colorIdx = 1 + round(alpha * (size(coolColormap, 1) - 1));
    color = coolColormap(colorIdx, :);

    t = markerT{i}(:);
    y = markerY{i}(:, 2);
    tRel = t - t(1);
    tStep = tRel(end);
    tLong = tRel;
    yLong = y;
    repeatCount = ceil((plotDuration + tStep) / tStep);
    for k = 1:repeatCount
        tLong = [tLong; tRel(2:end) + k * tStep]; %#ok<AGROW>
        yLong = [yLong; y(2:end)]; %#ok<AGROW>
    end
    i0 = find(islocalmax(yLong), 1, 'first');
    tLong = tLong(i0:end) - tLong(i0);
    yLong = yLong(i0:end);
    i1 = find(tLong >= plotDuration, 1, 'first');
    yEnd = interp1(tLong(i1 - 1:i1), yLong(i1 - 1:i1), plotDuration);
    tPlot = [tLong(1:i1 - 1); plotDuration];
    yPlot = [yLong(1:i1 - 1); yEnd];

    ax2 = nexttile(rightTiled, row);
    plot(ax2, tPlot, yPlot, 'LineWidth', 3, 'Color', color);
    grid(ax2, 'on');
    xlim(ax2, [0, plotDuration]);
    xticks(ax2, xTicks);
    ylim(ax2, concentrationLimits);
end

xlabel(ax2, 'Time (min)');
ylabel(rightTiled, 'concentration (a.u.)');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);
