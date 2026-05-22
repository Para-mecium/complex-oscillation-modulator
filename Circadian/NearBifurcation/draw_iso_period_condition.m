clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
circadianDir = fileparts(scriptDir);
repoDir = fileparts(circadianDir);
add_cbrewer_path(repoDir);
figureFile = fullfile(scriptDir, 'iso_period_condition.png');
curveDataDir = fullfile(scriptDir, 'data', 'iso_period', 'curves');
bifurcationFile = fullfile(scriptDir, 'data', 'bifurcation', 'circadian_bifurcation_line.mat');

%% Files
curvePeriods = [23.5, 24.0, 24.5 25];

curveFiles = cell(1, numel(curvePeriods));
for i = 1:numel(curvePeriods)
    curveFiles{i} = fullfile(curveDataDir, ...
        sprintf('iso_period_T%s_condition.mat', period_tag(curvePeriods(i))));
end

%% Plot settings
kdScale = 1e4;
xLimits = kdScale * [0, 1.5e-4];
yLimits = [0, 0.09];

nColor = 256;
conditionColorGamma = 0.85;
redLogCondition = 8.0;
colorbarMaxLogCondition = 10.0;
maxPlottedLogCondition = 9.8;
lineColor = [0.25, 0.25, 0.25];
lineStyle = '--';
lineWidth = 1;
maxLinePoints = 40;
labelOffsets = [ ...
    0.02, 0.000; ...
    0.02, -0.004; ...
    0.02, 0.000; ...
    0.02, 0.005];

%% Load bifurcation line
loadedBifurcation = load(bifurcationFile);
visibleHopfCurve = loadedBifurcation.visibleHopfCurve;
[hopfKdPlot, hopfATPlot] = downsample_line( ...
    visibleHopfCurve.Kd(:), visibleHopfCurve.AT(:), [], maxLinePoints);

%% Load iso-period curves
curveComponents = cell(1, numel(curveFiles));
curveAnchors = cell(1, numel(curveFiles));
allCurveLogConditions = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});

    if ~isfield(loaded, 'seed') || ~isfield(loaded, 'downBranch') || ~isfield(loaded, 'upBranch')
        error('draw_iso_period_condition:LegacyCurveCache', ...
            'Run build_iso_period_condition_data.m before drawing. Legacy four-branch cache is no longer supported: %s', ...
            curveFiles{i});
    end
    curveComponents{i} = join_at_branches(loaded.downBranch, loaded.seed, loaded.upBranch);
    curveAnchors{i} = curve_label_anchor_from_component(curveComponents{i});

    allCurveLogConditions = [allCurveLogConditions; curveComponents{i}.logCondition(:)]; %#ok<AGROW>
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
fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');
colormap(ax, conditionColormap);

for i = 1:numel(curveComponents)
    component = curveComponents{i};
    component = break_large_condition_points(component, maxPlottedLogCondition);
    component = break_out_of_range_points(component, kdScale, xLimits, yLimits);
    component = downsample_component(component, maxLinePoints);
    surface(ax, ...
        kdScale * [component.Kd.'; component.Kd.'], ...
        [component.AT.'; component.AT.'], ...
        zeros(2, numel(component.Kd)), ...
        [component.logCondition.'; component.logCondition.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1);

    anchor = curveAnchors{i};
    % text(ax, kdScale * anchor(1) + labelOffsets(i, 1), anchor(2) + labelOffsets(i, 2), ...
    %     sprintf('T=%.1f', curvePeriods(i)), ...
    %     'FontSize', 11, 'Color', [0.12, 0.32, 0.68]);
end

plot(ax, kdScale * hopfKdPlot, hopfATPlot, ...
    'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);

grid(ax, 'on');
box(ax, 'on');
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
% title(ax, 'Iso-period curves near bifurcation');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%%
function component = join_at_branches(downBranch, seed, upBranch)
downKd = flipud(branch_kd(downBranch));
downAT = flipud(branch_at(downBranch));
downLogCondition = flipud(branch_log_condition(downBranch));

upKd = branch_kd(upBranch);
upAT = branch_at(upBranch);
upLogCondition = branch_log_condition(upBranch);

component = struct();
component.Kd = [downKd; seed.Parameters(1); upKd];
component.AT = [downAT; seed.Parameters(2); upAT];
component.logCondition = fill_nonfinite_values([downLogCondition; log10(seed.conditionNumber); upLogCondition]);
end

function anchor = curve_label_anchor_from_component(component)
validMask = isfinite(component.Kd) & isfinite(component.AT);
if ~any(validMask)
    anchor = [0, 0];
    return
end

validKd = component.Kd(validMask);
validAT = component.AT(validMask);
[~, idx] = max(validKd);
anchor = [validKd(idx), validAT(idx)];
end

function values = branch_kd(branch)
if isempty(branch)
    values = zeros(0, 1);
    return
end
params = vertcat(branch.params);
values = params(:, 1);
end

function values = branch_at(branch)
if isempty(branch)
    values = zeros(0, 1);
    return
end
params = vertcat(branch.params);
values = params(:, 2);
end

function values = branch_log_condition(branch)
if isempty(branch)
    values = zeros(0, 1);
    return
end
values = log10(1 ./ [branch.directConditionEstimate]).';
values = fill_nonfinite_values(values);
end

function component = downsample_component(component, maxPoints)
[component.Kd, component.AT, component.logCondition] = downsample_line( ...
    component.Kd, component.AT, component.logCondition, maxPoints);
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

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end

function add_cbrewer_path(repoDir)
cbrewerDir = fullfile(repoDir, '+utils', 'cbrewer2');
if exist(cbrewerDir, 'dir') == 7 && ~contains(path, cbrewerDir)
    addpath(cbrewerDir);
end
end

function cmap = get_condition_colormap(numColors, gamma, colorRange, redStartPosition)
if exist('cbrewer2', 'file') ~= 2
    error('draw_iso_period_condition:MissingCbrewer2', ...
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

function component = break_out_of_range_points(component, kdScale, xLimits, yLimits)
plotX = kdScale * component.Kd(:);
plotY = component.AT(:);
inRange = plotX >= xLimits(1) & plotX <= xLimits(2) & ...
    plotY >= yLimits(1) & plotY <= yLimits(2);

newKd = zeros(0, 1);
newAT = zeros(0, 1);
newLogCondition = zeros(0, 1);

for i = 1:numel(component.Kd)
    if inRange(i)
        newKd = [newKd; component.Kd(i)]; %#ok<AGROW>
        newAT = [newAT; component.AT(i)]; %#ok<AGROW>
        newLogCondition = [newLogCondition; component.logCondition(i)]; %#ok<AGROW>
    elseif isempty(newKd) || ~isnan(newKd(end))
        newKd = [newKd; NaN]; %#ok<AGROW>
        newAT = [newAT; NaN]; %#ok<AGROW>
        newLogCondition = [newLogCondition; NaN]; %#ok<AGROW>
    end
end

component.Kd = newKd;
component.AT = newAT;
component.logCondition = newLogCondition;
end

function component = break_large_condition_points(component, maxLogCondition)
keepMask = component.logCondition(:) <= maxLogCondition;

newKd = zeros(0, 1);
newAT = zeros(0, 1);
newLogCondition = zeros(0, 1);

for i = 1:numel(component.Kd)
    if keepMask(i)
        newKd = [newKd; component.Kd(i)]; %#ok<AGROW>
        newAT = [newAT; component.AT(i)]; %#ok<AGROW>
        newLogCondition = [newLogCondition; component.logCondition(i)]; %#ok<AGROW>
    elseif isempty(newKd) || ~isnan(newKd(end))
        newKd = [newKd; NaN]; %#ok<AGROW>
        newAT = [newAT; NaN]; %#ok<AGROW>
        newLogCondition = [newLogCondition; NaN]; %#ok<AGROW>
    end
end

component.Kd = newKd;
component.AT = newAT;
component.logCondition = newLogCondition;
end

function values = fill_nonfinite_values(values)
finiteMask = isfinite(values);
if all(finiteMask)
    return
end
if ~any(finiteMask)
    values(:) = 0;
    return
end

idx = (1:numel(values)).';
finiteIdx = find(finiteMask);
values(1:finiteIdx(1) - 1) = values(finiteIdx(1));
values(finiteIdx(end) + 1:end) = values(finiteIdx(end));

finiteMask = isfinite(values);
values(~finiteMask) = interp1(idx(finiteMask), values(finiteMask), idx(~finiteMask), 'linear');
end
