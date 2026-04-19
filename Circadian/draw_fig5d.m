clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
curveDataDir = fullfile(scriptDir, 'data', 'fig5d', 'curves');
markerDataDir = fullfile(scriptDir, 'data', 'fig5d', 'markers');
figureFile = fullfile(scriptDir, 'fig5d.png');

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

for i = 1:numel(markerFiles)
    loaded = load(markerFiles{i}, 'Parameters', 'period');
    markerParams(i, :) = loaded.Parameters(:, 1:2);
    markerPeriodValues(i) = loaded.period;
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 640, 520]);
ax = axes(fig);
hold(ax, 'on');
colormap(ax, coolColormap);
for i = 1:numel(curveData)
    curve = curveData{i};
    surface(ax, ...
        kdScale * [curve.Kd.'; curve.Kd.'], ...
        [curve.AT.'; curve.AT.'], ...
        zeros(2, numel(curve.Kd)), ...
        [curve.period.'; curve.period.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);
end
scatter(ax, kdScale * markerParams(:, 1), markerParams(:, 2), 120, markerPeriodValues, ...
    'filled', 'MarkerEdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);
grid(ax, 'on');
clim(ax, periodLimits);
cb = colorbar(ax);
cb.Label.String = 'Period (hour)';
xlabel(ax, 'K_d (\times 10^{-4})');
ylabel(ax, 'A_T (a.u.)');
xlim(ax, xLimits);
ylim(ax, yLimits);
title(ax, 'Fig. 5d: Iso-maximum curves');

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

function tag = maximum_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end
