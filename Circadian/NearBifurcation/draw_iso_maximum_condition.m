clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
circadianDir = fileparts(scriptDir);
repoDir = fileparts(circadianDir);
add_cbrewer_path(repoDir);
figureFile = fullfile(scriptDir, 'iso_maximum_condition.png');
curveDataDir = fullfile(scriptDir, 'data', 'iso_maximum', 'curves');
bifurcationFile = fullfile(scriptDir, 'data', 'bifurcation', 'circadian_bifurcation_line.mat');

%% Files
targetMaxima = [0.20, 0.16, 0.12, 0.08, 0.04];

curveFiles = cell(1, numel(targetMaxima));
for i = 1:numel(targetMaxima)
    curveFiles{i} = fullfile(curveDataDir, ...
        sprintf('iso_maximum_M%s_condition.mat', maximum_tag(targetMaxima(i))));
end

%% Plot settings
kdScale = 1e4;
xLimits = kdScale * [0, 1.5e-4];
yLimits = [0, 0.09];

nColor = 256;
conditionColorGamma = 0.85;
redLogCondition = 8.0;
colorbarMaxLogCondition = 10.0;
lineColor = [0.25, 0.25, 0.25];
lineStyle = '--';
lineWidth = 1.4;

%% Load bifurcation line
loadedBifurcation = load(bifurcationFile);
visibleHopfCurve = loadedBifurcation.visibleHopfCurve;

%% Load iso-maximum curves
curveData = cell(1, numel(curveFiles));
allCurveLogConditions = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});
    curveData{i} = join_curve_branches(loaded.leftBranch, loaded.seed, loaded.rightBranch);
    allCurveLogConditions = [allCurveLogConditions; curveData{i}.logCondition(:)]; %#ok<AGROW>
end

conditionLimits = [min(allCurveLogConditions), max(allCurveLogConditions)];
if conditionLimits(1) >= conditionLimits(2)
    delta = max(abs(conditionLimits(1)) * 0.05, 1e-6);
    conditionLimits = [conditionLimits(1) - delta, conditionLimits(2) + delta];
end
colorConditionLimits = [conditionLimits(1), colorbarMaxLogCondition];
if colorConditionLimits(1) >= colorConditionLimits(2)
    colorConditionLimits(2) = colorConditionLimits(1) + 1e-6;
end
redStartPosition = (redLogCondition - colorConditionLimits(1)) / diff(colorConditionLimits);
conditionColormap = get_condition_colormap(nColor, conditionColorGamma, [0.04, 0.96], redStartPosition);

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 640, 520]);
ax = axes(fig);
hold(ax, 'on');
colormap(ax, conditionColormap);
for i = 1:numel(curveData)
    curve = curveData{i};
    surface(ax, ...
        kdScale * [curve.Kd.'; curve.Kd.'], ...
        [curve.AT.'; curve.AT.'], ...
        zeros(2, numel(curve.Kd)), ...
        [curve.logCondition.'; curve.logCondition.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
end
plot(ax, kdScale * visibleHopfCurve.Kd(:), visibleHopfCurve.AT(:), ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
grid(ax, 'on');
clim(ax, colorConditionLimits);
cb = colorbar(ax);
cb.Label.String = 'log_{10} condition number';
conditionTicks = condition_colorbar_ticks(colorConditionLimits);
conditionTicks = conditionTicks(conditionTicks <= colorbarMaxLogCondition);
cb.Ticks = conditionTicks(:);
cb.TickLabels = arrayfun(@(value) sprintf('%.1f', value), conditionTicks(:), 'UniformOutput', false);
xlabel(ax, 'K_d (\times 10^{-4})');
ylabel(ax, 'A_T (a.u.)');
xlim(ax, xLimits);
ylim(ax, yLimits);
title(ax, 'Iso-maximum curves near bifurcation');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%%
function curve = join_curve_branches(leftBranch, seed, rightBranch)
if isempty(leftBranch)
    leftKd = zeros(0, 1);
    leftAT = zeros(0, 1);
    leftLogCondition = zeros(0, 1);
else
    leftParameters = vertcat(leftBranch.params);
    leftKd = flipud(leftParameters(:, 1));
    leftAT = flipud(leftParameters(:, 2));
    leftLogCondition = flipud(branch_log_condition(leftBranch));
end

if isempty(rightBranch)
    rightKd = zeros(0, 1);
    rightAT = zeros(0, 1);
    rightLogCondition = zeros(0, 1);
else
    rightParameters = vertcat(rightBranch.params);
    rightKd = rightParameters(:, 1);
    rightAT = rightParameters(:, 2);
    rightLogCondition = branch_log_condition(rightBranch);
end

curve = struct();
curve.Kd = [leftKd; seed.Parameters(1); rightKd];
curve.AT = [leftAT; seed.Parameters(2); rightAT];
curve.logCondition = [leftLogCondition; log10(seed.conditionNumber); rightLogCondition];
end

function values = branch_log_condition(branch)
values = log10(1 ./ [branch.directConditionEstimate]).';
end

function tag = maximum_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function add_cbrewer_path(repoDir)
cbrewerDir = fullfile(repoDir, '+utils', 'cbrewer2');
if exist(cbrewerDir, 'dir') == 7 && ~contains(path, cbrewerDir)
    addpath(cbrewerDir);
end
end

function cmap = get_condition_colormap(numColors, gamma, colorRange, redStartPosition)
if exist('cbrewer2', 'file') ~= 2
    error('draw_iso_maximum_condition:MissingCbrewer2', ...
        'Could not find cbrewer2. Expected +utils/cbrewer2/cbrewer2.m under the repository root.');
end

if nargin < 3
    colorRange = [0, 1];
end
if nargin < 4
    redStartPosition = 1;
end

basePosition = linspace(0, 1, numColors).';
baseMap = flipud(cbrewer2('RdYlBu', numColors, 'linear', 'rgb'));
redStartPosition = min(max(redStartPosition, eps), 1);
redColorPosition = 0.94;
warpedPosition = zeros(size(basePosition));
belowMask = basePosition <= redStartPosition;
warpedPosition(belowMask) = redColorPosition * (basePosition(belowMask) / redStartPosition) .^ gamma;
if redStartPosition < 1
    upperPosition = (basePosition(~belowMask) - redStartPosition) / (1 - redStartPosition);
    warpedPosition(~belowMask) = redColorPosition + (1 - redColorPosition) * upperPosition;
else
    warpedPosition(~belowMask) = 1;
end
warpedPosition = colorRange(1) + diff(colorRange) * warpedPosition;
cmap = interp1(basePosition, baseMap, warpedPosition, 'linear');

middleWeight = exp(-((basePosition - 0.5) / 0.18) .^ 2);
middleBlendStrength = 0.60;
middleColor = [0.88, 0.58, 0.18];
cmap = (1 - middleBlendStrength * middleWeight) .* cmap + ...
    (middleBlendStrength * middleWeight) .* middleColor;
cmap = min(max(cmap, 0), 1);
end

function ticks = condition_colorbar_ticks(limits)
tickStep = 0.5;
firstTick = ceil(limits(1) / tickStep) * tickStep;
lastTick = floor(limits(2) / tickStep) * tickStep;
ticks = firstTick:tickStep:lastTick;
if isempty(ticks)
    ticks = linspace(limits(1), limits(2), 5);
end
end
