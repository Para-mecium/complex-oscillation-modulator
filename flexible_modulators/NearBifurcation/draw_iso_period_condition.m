clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
flexmodDir = fileparts(scriptDir);
figureFile = fullfile(scriptDir, 'iso_period_condition.png');
bifurcationFile = fullfile(flexmodDir, 'data', 'fig3', 'fig3_bifurcation_line.mat');
curveDataDir = fullfile(scriptDir, 'data', 'iso_period', 'curves');

%% Files
curvePeriods = [50, 60, 70, 80, 90];

curveFiles = cell(1, numel(curvePeriods));
for i = 1:numel(curvePeriods)
    curveFiles{i} = fullfile(curveDataDir, sprintf('iso_period_T%03d_condition.mat', round(curvePeriods(i))));
end

%% Plot settings
lineColor = [0.25, 0.25, 0.25];
lineStyle = '--';
lineWidth = 1.4;

nColor = 256;
conditionColormap = parula(nColor);

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
curveLogCondition = cell(1, numel(curveFiles));
allCurveLogConditions = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});

    seedParams = loaded.seed.params(:).';
    seedLogCondition = log10(loaded.seed.conditionNumber);

    leftParams = vertcat(loaded.leftBranch.params);
    leftCondition = 1 ./ [loaded.leftBranch.directConditionEstimate].';
    rightParams = vertcat(loaded.rightBranch.params);
    rightCondition = 1 ./ [loaded.rightBranch.directConditionEstimate].';

    curveI1{i} = [flipud(leftParams(:, 1)); seedParams(1); rightParams(:, 1)];
    curveET{i} = [flipud(leftParams(:, 2)); seedParams(2); rightParams(:, 2)];
    curveLogCondition{i} = [flipud(log10(leftCondition)); seedLogCondition; log10(rightCondition)];

    allCurveLogConditions = [allCurveLogConditions; curveLogCondition{i}]; %#ok<AGROW>
end

conditionLimits = [min(allCurveLogConditions), max(allCurveLogConditions)];
if conditionLimits(1) >= conditionLimits(2)
    delta = max(abs(conditionLimits(1)) * 0.05, 1e-6);
    conditionLimits = [conditionLimits(1) - delta, conditionLimits(2) + delta];
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 720, 560]);
ax = axes(fig);
hold(ax, 'on');
colormap(ax, conditionColormap);

for i = 1:numel(curveFiles)
    surface(ax, [curveI1{i}.'; curveI1{i}.'], [curveET{i}.'; curveET{i}.'], ...
        zeros(2, numel(curveI1{i})), [curveLogCondition{i}.'; curveLogCondition{i}.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);

    text(ax, curveI1{i}(end) - 0.02, curveET{i}(end), sprintf('%d', curvePeriods(i)), ...
        'FontSize', 11, 'Color', [0.15, 0.35, 0.75], 'HorizontalAlignment', 'right');
end

plot(ax, hopfI1Plot, hopfETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
plot(ax, lpI1Plot, lpETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);

grid(ax, 'on');
clim(ax, conditionLimits);
cb = colorbar(ax);
cb.Label.String = 'log_{10} condition number';
xlabel(ax, 'I_1 (a.u.)');
ylabel(ax, 'E_T (a.u.)');
title(ax, 'Iso-period curves near bifurcation');
xlim(ax, xLimits);
ylim(ax, yLimits);

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);
