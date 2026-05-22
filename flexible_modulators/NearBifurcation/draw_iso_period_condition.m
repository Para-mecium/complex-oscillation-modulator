clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
flexmodDir = fileparts(scriptDir);
repoDir = fileparts(flexmodDir);
add_cbrewer_path(repoDir);
figureFile = fullfile(scriptDir, 'iso_period_condition.png');
bifurcationFile = fullfile(flexmodDir, 'data', 'fig3', 'fig3_bifurcation_line.mat');
curveDataDir = fullfile(scriptDir, 'data', 'iso_period', 'curves');

%% Files
curvePeriods = [50, 60, 70, 80, 90];
referenceAmplitudes = 1.2:0.3:3.6;
referenceCurveDataDir = fullfile(scriptDir, 'data', 'iso_amplitude', 'curves');

curveFiles = cell(1, numel(curvePeriods));
for i = 1:numel(curvePeriods)
    curveFiles{i} = fullfile(curveDataDir, sprintf('iso_period_T%03d_condition.mat', round(curvePeriods(i))));
end

referenceCurveFiles = cell(1, numel(referenceAmplitudes));
for i = 1:numel(referenceAmplitudes)
    amplitudeTag = strrep(sprintf('%.1f', referenceAmplitudes(i)), '.', 'p');
    referenceCurveFiles{i} = fullfile(referenceCurveDataDir, sprintf('iso_amplitude_A%s_condition.mat', amplitudeTag));
end

%% Plot settings
lineColor = [0.25, 0.25, 0.25];
lineStyle = '--';
lineWidth = 1;
maxLinePoints = 40;

nColor = 256;
% gamma < 1 expands low-to-mid condition values and compresses the high-end tail.
conditionColorGamma = 0.45;
conditionColormap = get_condition_colormap(nColor, conditionColorGamma);

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
[hopfI1Plot, hopfETPlot] = downsample_line(hopfI1Plot, hopfETPlot, [], maxLinePoints);

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
[lpI1Plot, lpETPlot] = downsample_line(lpI1Plot, lpETPlot, [], maxLinePoints);

%% Load iso-period curves
curveI1 = cell(1, numel(curveFiles));
curveET = cell(1, numel(curveFiles));
curveLogCondition = cell(1, numel(curveFiles));

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
    [curveI1{i}, curveET{i}, curveLogCondition{i}] = downsample_line( ...
        curveI1{i}, curveET{i}, curveLogCondition{i}, maxLinePoints);
end

conditionLimits = get_curve_log_condition_limits(referenceCurveFiles);

%% Draw figure
fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');
colormap(ax, conditionColormap);

for i = 1:numel(curveFiles)
    surface(ax, [curveI1{i}.'; curveI1{i}.'], [curveET{i}.'; curveET{i}.'], ...
        zeros(2, numel(curveI1{i})), [curveLogCondition{i}.'; curveLogCondition{i}.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1);

    % text(ax, curveI1{i}(end) - 0.02, curveET{i}(end), sprintf('%d', curvePeriods(i)), ...
    %     'FontSize', 11, 'Color', [0.15, 0.35, 0.75], 'HorizontalAlignment', 'right');
end

plot(ax, hopfI1Plot, hopfETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
plot(ax, lpI1Plot, lpETPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);

grid(ax, 'on');
box(ax, 'on');
clim(ax, conditionLimits);
cb = colorbar(ax);
cb.Label.String = 'log_{10} condition number';
xlabel(ax, 'I_1 (a.u.)');
ylabel(ax, 'E_T (a.u.)');
% title(ax, 'Iso-period curves near bifurcation');
xlim(ax, xLimits);
ylim(ax, yLimits);

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

function add_cbrewer_path(repoDir)
cbrewerDir = fullfile(repoDir, '+utils', 'cbrewer2');
if exist(cbrewerDir, 'dir') == 7 && ~contains(path, cbrewerDir)
    addpath(cbrewerDir);
end
end

function cmap = get_condition_colormap(numColors, gamma)
if nargin < 1
    numColors = 256;
end
if nargin < 2
    gamma = 1;
end
if exist('cbrewer2', 'file') ~= 2
    error('draw_iso_period_condition:MissingCbrewer2', ...
        'Could not find cbrewer2. Expected +utils/cbrewer2/cbrewer2.m under the repository root.');
end
baseMap = flipud(cbrewer2('RdYlBu', numColors, 'linear', 'rgb'));
basePosition = linspace(0, 1, numColors).';
warpedPosition = basePosition .^ gamma;
cmap = interp1(basePosition, baseMap, warpedPosition, 'linear');
end

function conditionLimits = get_curve_log_condition_limits(curveFiles)
allCurveLogConditions = [];
for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});

    seedLogCondition = log10(loaded.seed.conditionNumber);
    leftCondition = 1 ./ [loaded.leftBranch.directConditionEstimate].';
    rightCondition = 1 ./ [loaded.rightBranch.directConditionEstimate].';
    curveLogCondition = [log10(leftCondition); seedLogCondition; log10(rightCondition)];

    allCurveLogConditions = [allCurveLogConditions; curveLogCondition]; %#ok<AGROW>
end

conditionLimits = [min(allCurveLogConditions), max(allCurveLogConditions)];
if conditionLimits(1) >= conditionLimits(2)
    delta = max(abs(conditionLimits(1)) * 0.05, 1e-6);
    conditionLimits = [conditionLimits(1) - delta, conditionLimits(2) + delta];
end
end

function [xOut, yOut, cOut] = downsample_line(x, y, c, maxPoints)
x = x(:);
y = y(:);
hasColor = ~isempty(c);
if hasColor
    c = c(:);
    finiteMask = isfinite(x) & isfinite(y) & isfinite(c);
else
    finiteMask = isfinite(x) & isfinite(y);
end
finiteIdx = find(finiteMask);
if isempty(finiteIdx)
    xOut = zeros(0, 1);
    yOut = zeros(0, 1);
    cOut = c;
    return
end
if numel(finiteIdx) == 1
    xOut = repmat(x(finiteIdx), maxPoints, 1);
    yOut = repmat(y(finiteIdx), maxPoints, 1);
    if hasColor
        cOut = repmat(c(finiteIdx), maxPoints, 1);
    else
        cOut = [];
    end
    return
end
if numel(finiteIdx) < maxPoints
    samplePosition = linspace(1, numel(finiteIdx), maxPoints).';
    finitePosition = (1:numel(finiteIdx)).';
    xOut = interp1(finitePosition, x(finiteIdx), samplePosition, 'linear');
    yOut = interp1(finitePosition, y(finiteIdx), samplePosition, 'linear');
    if hasColor
        cOut = interp1(finitePosition, c(finiteIdx), samplePosition, 'linear');
    else
        cOut = [];
    end
    return
end
keepFiniteIdx = finiteIdx(round(linspace(1, numel(finiteIdx), maxPoints)));
keepMask = false(size(x));
keepMask(keepFiniteIdx) = true;
nanIdx = find(~finiteMask);
for i = 1:numel(nanIdx)
    hasKeptBefore = any(keepFiniteIdx < nanIdx(i));
    hasKeptAfter = any(keepFiniteIdx > nanIdx(i));
    if hasKeptBefore && hasKeptAfter
        keepMask(nanIdx(i)) = true;
    end
end
xOut = x(keepMask);
yOut = y(keepMask);
if hasColor
    cOut = c(keepMask);
else
    cOut = [];
end
end
