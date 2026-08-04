clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
curveDataDir = fullfile(scriptDir, 'data', 'fig5d', 'curves');
markerDataDir = fullfile(scriptDir, 'data', 'fig5d', 'markers');
figureFile = fullfile(scriptDir, 'figS15b.png');

%% Files
targetMaxima = [0.20, 0.16, 0.12, 0.08, 0.04];
markerMaximum = 0.12;
markerPeriods = [23.5, 24.0, 24.5];

curveFiles = cell(1, numel(targetMaxima));
for i = 1:numel(targetMaxima)
    curveFiles{i} = fullfile(curveDataDir, ...
        sprintf('fig5d_iso_maximum_curve_M%s.mat', maximum_tag(targetMaxima(i))));
end

markerFiles = cell(1, numel(markerPeriods));
for i = 1:numel(markerPeriods)
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('fig5d_marker_M%s_T%s.mat', ...
        maximum_tag(markerMaximum), period_tag(markerPeriods(i))));
end

%% Plot settings
kdScale = 1e4;
xLimits = kdScale * [0, 1.5e-4];
yLimits = [0.01, 0.09];

nColor = 256;
tColor = linspace(0, 1, nColor).';
coolColormap = (1 - tColor) .* [0.05, 0.25, 0.60] + tColor .* [0.92, 0.95, 0.98];

%% Load iso-maximum curves
curveData = cell(1, numel(curveFiles));
allCurvePeriods = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});
    curveData{i} = join_curve_branches(loaded.leftBranch, loaded.seed, loaded.rightBranch);
    allCurvePeriods = [allCurvePeriods; curveData{i}.period(:)]; %#ok<AGROW>
end

periodLimits = [min(allCurvePeriods), max(allCurvePeriods)];

%% Load marker data
markerParams = zeros(numel(markerFiles), 2);
markerPeriodValues = zeros(numel(markerFiles), 1);
shiftedTSeries = cell(1, numel(markerFiles));
shiftedObsSeries = cell(1, numel(markerFiles));

for i = 1:numel(markerFiles)
    loaded = load(markerFiles{i}, 'Parameters', 'TS', 'period');
    markerParams(i, :) = loaded.Parameters(:, 1:2);
    markerPeriodValues(i) = loaded.period;

    t = loaded.TS{1};
    y = loaded.TS{2};
    obs = y(:, 2) + y(:, 3);
    [shiftedTSeries{i}, shiftedObsSeries{i}] = shift_cycle_to_max_local(t, y, obs, loaded.period(1));
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 1100, 460]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on');
colormap(ax1, coolColormap);

for i = 1:numel(curveData)
    curve = curveData{i};
    surface(ax1, ...
        kdScale * [curve.Kd.'; curve.Kd.'], ...
        [curve.AT.'; curve.AT.'], ...
        zeros(2, numel(curve.Kd)), ...
        [curve.period.'; curve.period.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
end

scatter(ax1, kdScale * markerParams(:, 1), markerParams(:, 2), 120, markerPeriodValues, ...
    'filled', 'MarkerEdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);

grid(ax1, 'on');
clim(ax1, periodLimits);
cb = colorbar(ax1);
cb.Label.String = 'Period (hour)';
xlabel(ax1, 'K_d (\times 10^{-4})');
ylabel(ax1, 'A_T (a.u.)');
xlim(ax1, xLimits);
ylim(ax1, yLimits);
title(ax1, 'Fig. S15b1: Iso-maximum curves');

ax2 = nexttile(tl, 2);
hold(ax2, 'on');

for i = 1:numel(markerFiles)
    lineColor = value_to_color(markerPeriodValues(i), periodLimits, coolColormap);
    plot(ax2, shiftedTSeries{i}, shiftedObsSeries{i}, 'LineWidth', 2.5, ...
        'Color', lineColor, 'DisplayName', sprintf('T = %.1f', markerPeriods(i)));
end

grid(ax2, 'on');
xlabel(ax2, 'Time (hour)');
ylabel(ax2, 'P_{tot} (a.u.)');
title(ax2, sprintf('Fig. S15b2: Time series (max = %.3f)', markerMaximum));
legend(ax2, 'Location', 'best');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%%
function curve = join_curve_branches(leftBranch, seed, rightBranch)
if isempty(leftBranch)
    leftKd = zeros(0, 1);
    leftAT = zeros(0, 1);
    leftPeriod = zeros(0, 1);
else
    leftParameters = vertcat(leftBranch.params);
    leftKd = flipud(leftParameters(:, 1));
    leftAT = flipud(leftParameters(:, 2));
    leftPeriod = flipud(branch_period(leftBranch));
end

if isempty(rightBranch)
    rightKd = zeros(0, 1);
    rightAT = zeros(0, 1);
    rightPeriod = zeros(0, 1);
else
    rightParameters = vertcat(rightBranch.params);
    rightKd = rightParameters(:, 1);
    rightAT = rightParameters(:, 2);
    rightPeriod = branch_period(rightBranch);
end

curve = struct();
curve.Kd = [leftKd; seed.Parameters(1); rightKd];
curve.AT = [leftAT; seed.Parameters(2); rightAT];
curve.period = [leftPeriod; seed.period; rightPeriod];
end

function values = branch_period(branch)
values = arrayfun(@(entry) entry.derived.period, branch).';
end

function [shiftedT, shiftedObs] = shift_cycle_to_max_local(t, y, obs, period)
obs = obs(:) + 0 .* (y(:, 2) + y(:, 3));
[~, peakIdx] = max(obs);
shiftedT = mod(t(:) - t(peakIdx), period);
[shiftedT, order] = sort(shiftedT);
shiftedObs = obs(order);
end

function color = value_to_color(value, valueLimits, cmap)
position = 1 + (size(cmap, 1) - 1) * (value - valueLimits(1)) / (valueLimits(2) - valueLimits(1));
lowerIdx = floor(position);
upperIdx = ceil(position);
weight = position - lowerIdx;
color = (1 - weight) * cmap(lowerIdx, :) + weight * cmap(upperIdx, :);
end

function tag = maximum_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end
