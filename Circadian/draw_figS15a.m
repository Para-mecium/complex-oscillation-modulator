clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
figureFile = fullfile(scriptDir, 'figS15a.png');
curveDataDir = fullfile(scriptDir, 'data', 'figS15a', 'curves');
markerDataDir = fullfile(scriptDir, 'data', 'figS15a', 'markers');

%% Files
targetAmplitudes = [0.005, 0.015, 0.025, 0.035, 0.0385];
markerAmplitude = 0.025;
markerPeriods = [23.5, 24.0, 24.5];

curveFiles = cell(1, numel(targetAmplitudes));
for i = 1:numel(targetAmplitudes)
    curveFiles{i} = fullfile(curveDataDir, ...
        sprintf('figS15a_iso_amplitude_curve_A%s.mat', amplitude_tag(targetAmplitudes(i))));
end

markerFiles = cell(1, numel(markerPeriods));
for i = 1:numel(markerPeriods)
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('figS15a_marker_A%s_T%s.mat', amplitude_tag(markerAmplitude), period_tag(markerPeriods(i))));
end

%% Plot settings
kdScale = 1e4;
xLimits = kdScale * [0, 1.5e-4];
yLimits = [0.01, 0.09];

nColor = 256;
tColor = linspace(0, 1, nColor).';
coolColormap = (1 - tColor) .* [0.05, 0.25, 0.60] + tColor .* [0.92, 0.95, 0.98];

%% Load iso-amplitude curves
curveData = cell(1, numel(curveFiles));
allCurvePeriods = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i}, 'seed', 'leftBranch', 'rightBranch');
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

    markerParams(i, :) = loaded.Parameters(1, 1:2);
    markerPeriodValues(i) = loaded.period(1);

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
    component = curveData{i};
    surface(ax1, ...
        kdScale * [component.Kd.'; component.Kd.'], ...
        [component.AT.'; component.AT.'], ...
        zeros(2, numel(component.Kd)), ...
        [component.period.'; component.period.'], ...
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
title(ax1, 'Fig. S15a1: Iso-amplitude curves');

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
title(ax2, sprintf('Fig. S15a2: Time series (A = %.3f)', markerAmplitude));
legend(ax2, 'Location', 'best');

exportgraphics(fig, figureFile, 'Resolution', 300);

%%
function component = join_curve_branches(leftBranch, seed, rightBranch)
if isempty(leftBranch)
    leftParameters = zeros(0, numel(seed.Parameters));
    leftPeriods = zeros(0, 1);
else
    leftParameters = vertcat(leftBranch.params);
    leftPeriods = branch_period(leftBranch);
end

if isempty(rightBranch)
    rightParameters = zeros(0, numel(seed.Parameters));
    rightPeriods = zeros(0, 1);
else
    rightParameters = vertcat(rightBranch.params);
    rightPeriods = branch_period(rightBranch);
end

component = struct();
component.Kd = [flipud(leftParameters(:, 1)); seed.Parameters(1); rightParameters(:, 1)];
component.AT = [flipud(leftParameters(:, 2)); seed.Parameters(2); rightParameters(:, 2)];
component.period = [flipud(leftPeriods(:)); seed.period(1); rightPeriods(:)];
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

function tag = amplitude_tag(value)
tag = regexprep(sprintf('%.4f', value), '0+$', '');
tag = regexprep(tag, '\.$', '');
tag = strrep(tag, '.', 'p');
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end
