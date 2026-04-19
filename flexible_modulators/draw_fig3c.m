clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
figureFile = fullfile(scriptDir, 'fig3c.png');
bifurcationFile = fullfile(scriptDir, 'data', 'fig3', 'fig3_bifurcation_line.mat');
curveDataDir = fullfile(scriptDir, 'data', 'fig3c', 'curves');
markerDataDir = fullfile(scriptDir, 'data', 'fig3c', 'markers');

%% Files
curvePeriods = [50, 60, 70, 80, 90];
markerPeriod = 80;
markerAmplitudes = [2.0, 2.5, 3.0];

curveFiles = cell(1, numel(curvePeriods));
for i = 1:numel(curvePeriods)
    curveFiles{i} = fullfile(curveDataDir, sprintf('fig3c_iso_period_T%03d.mat', round(curvePeriods(i))));
end

markerFiles = cell(1, numel(markerAmplitudes));
for i = 1:numel(markerAmplitudes)
    amplitudeTag = strrep(sprintf('%.1f', markerAmplitudes(i)), '.', 'p');
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('fig3c_marker_T%03d_A%s.mat', round(markerPeriod), amplitudeTag));
end

%% Plot settings
lineColor = [0.25, 0.25, 0.25];
lineStyle = '--';
lineWidth = 1.4;

nColor = 256;
tColor = linspace(0, 1, nColor).';
warmColormap = (1 - tColor) .* [1.0, 0.65, 0.22] + tColor .* [0.98, 0.96, 0.92];

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
%% Load iso-period curves
curveI1 = cell(1, numel(curveFiles));
curveET = cell(1, numel(curveFiles));
curveAmplitude = cell(1, numel(curveFiles));
allCurveAmplitudes = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});

    seedParams = loaded.seed.params(:).';
    seedAmplitude = loaded.seed.derivedView.varAmp(2);

    leftParams = vertcat(loaded.leftBranch.params);
    leftAmplitude = arrayfun(@(entry) entry.derived.varAmp(2), loaded.leftBranch).';
    rightParams = vertcat(loaded.rightBranch.params);
    rightAmplitude = arrayfun(@(entry) entry.derived.varAmp(2), loaded.rightBranch).';

    curveI1{i} = [flipud(leftParams(:, 1)); seedParams(1); rightParams(:, 1)];
    curveET{i} = [flipud(leftParams(:, 2)); seedParams(2); rightParams(:, 2)];
    curveAmplitude{i} = [flipud(leftAmplitude); seedAmplitude; rightAmplitude];

    allCurveAmplitudes = [allCurveAmplitudes; curveAmplitude{i}]; %#ok<AGROW>
end

ampLimits = [min(allCurveAmplitudes), max(allCurveAmplitudes)];
if ampLimits(1) >= ampLimits(2)
    delta = max(abs(ampLimits(1)) * 0.05, 1e-6);
    ampLimits = [ampLimits(1) - delta, ampLimits(2) + delta];
end

%% Load marker data
markerParams = cell(1, numel(markerFiles));
markerT = cell(1, numel(markerFiles));
markerY = cell(1, numel(markerFiles));
markerAmplitude = zeros(1, numel(markerFiles));

for i = 1:numel(markerFiles)
    loaded = load(markerFiles{i});
    markerParams{i} = loaded.Parameters(:).';
    markerT{i} = loaded.TS{1};
    markerY{i} = loaded.TS{2};
    markerAmplitude(i) = loaded.varAmp(2);
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 480]);
tiled = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled, 1);
hold(ax1, 'on');
colormap(ax1, warmColormap);

for i = 1:numel(curveFiles)
    surface(ax1, [curveI1{i}.'; curveI1{i}.'], [curveET{i}.'; curveET{i}.'], ...
        zeros(2, numel(curveI1{i})), [curveAmplitude{i}.'; curveAmplitude{i}.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);

    text(ax1, curveI1{i}(end) + 0.02, curveET{i}(end), sprintf('%d', curvePeriods(i)), ...
        'FontSize', 11, 'Color', [0.15, 0.35, 0.75]);
end

for i = 1:numel(markerFiles)
    scatter(ax1, markerParams{i}(1), markerParams{i}(2), 160, markerAmplitude(i), 'filled', ...
        'MarkerEdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 1.2);
end

plot(ax1, hopfI1Plot, hopfETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
plot(ax1, lpI1Plot, lpETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);

grid(ax1, 'on');
clim(ax1, ampLimits);
cb = colorbar(ax1);
cb.Label.String = 'Protein amplitude';
xlabel(ax1, 'I_1 (a.u.)');
ylabel(ax1, 'E_T (a.u.)');
title(ax1, sprintf('Iso-period curves (T = %d highlighted)', markerPeriod));
xlim(ax1, xLimits);
ylim(ax1, yLimits);

ax2 = nexttile(tiled, 2);
hold(ax2, 'on');
colormap(ax2, warmColormap);

for i = 1:numel(markerFiles)
    alpha = (markerAmplitude(i) - ampLimits(1)) / (ampLimits(2) - ampLimits(1));
    alpha = min(max(alpha, 0), 1);
    colorIdx = 1 + round(alpha * (size(warmColormap, 1) - 1));
    color = warmColormap(colorIdx, :);

    t = markerT{i}(:);
    y = markerY{i};
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

    plot(ax2, tPlot, yPlot, 'LineWidth', 3, 'Color', color, ...
        'DisplayName', sprintf('A = %.1f', markerAmplitudes(i)));
end

grid(ax2, 'on');
xlabel(ax2, 'Time (min)');
ylabel(ax2, 'Y (a.u.)');
title(ax2, sprintf('Period = %d min', markerPeriod));
legend(ax2, 'Location', 'northwest');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);
