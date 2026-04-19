clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
figureFile = fullfile(scriptDir, 'fig5b.png');
curveDataDir = fullfile(scriptDir, 'data', 'fig5b', 'curves');
markerDataDir = fullfile(scriptDir, 'data', 'fig5b', 'markers');

%% Files
curvePeriods = [23.5, 24.0, 24.5, 25.0];
markerPeriod = 24.0;
markerAmplitudes = [0.01, 0.02, 0.04];

curveFiles = cell(1, numel(curvePeriods));
for i = 1:numel(curvePeriods)
    curveFiles{i} = fullfile(curveDataDir, ...
        sprintf('fig5b_iso_period_curve_T%s.mat', period_tag(curvePeriods(i))));
end

markerFiles = cell(1, numel(markerAmplitudes));
for i = 1:numel(markerAmplitudes)
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('fig5b_marker_T%s_A%s.mat', period_tag(markerPeriod), amplitude_tag(markerAmplitudes(i))));
end

%% Plot settings
kdScale = 1e4;
xLimits = kdScale * [0, 1.5e-4];
yLimits = [0.01, 0.09];

nColor = 256;
tColor = linspace(0, 1, nColor).';
warmColormap = (1 - tColor) .* [1.00, 0.55, 0.14] + tColor .* [0.98, 0.96, 0.92];

%% Load iso-period curves
curveComponents = cell(1, numel(curveFiles));
curveAnchors = cell(1, numel(curveFiles));
allCurveAmplitudes = [];

for i = 1:numel(curveFiles)
    loaded = load(curveFiles{i});

    lower = join_curve_branches(loaded.lowerLeftBranch, loaded.lowerSeed, loaded.lowerRightBranch);
    upper = join_curve_branches(loaded.upperLeftBranch, loaded.upperSeed, loaded.upperRightBranch);
    curveComponents{i} = stitch_nearest_components(lower, upper);
    curveAnchors{i} = curve_label_anchor(loaded);

    allCurveAmplitudes = [allCurveAmplitudes; curveComponents{i}.obsAmp(:)]; %#ok<AGROW>
end

ampLimits = [min(allCurveAmplitudes), max(allCurveAmplitudes)];

%% Load marker data
markerParams = zeros(numel(markerFiles), 2);
markerAmplitudeValues = zeros(numel(markerFiles), 1);

for i = 1:numel(markerFiles)
    loaded = load(markerFiles{i});
    markerParams(i, :) = loaded.Parameters(:).';
    markerAmplitudeValues(i) = loaded.obsAmp(1);
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 640, 520]);
ax = axes(fig);
hold(ax, 'on');
colormap(ax, warmColormap);

for i = 1:numel(curveComponents)
    component = curveComponents{i};
    surface(ax, ...
        kdScale * [component.Kd.'; component.Kd.'], ...
        [component.AT.'; component.AT.'], ...
        zeros(2, numel(component.Kd)), ...
        [component.obsAmp.'; component.obsAmp.'], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 3);

    anchor = curveAnchors{i};
    text(ax, kdScale * anchor(1) + 0.02, anchor(2), sprintf('T=%.1f', curvePeriods(i)), ...
        'FontSize', 11, 'Color', [0.12, 0.32, 0.68]);
end

scatter(ax, kdScale * markerParams(:, 1), markerParams(:, 2), 120, markerAmplitudeValues, ...
    'filled', 'MarkerEdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);

grid(ax, 'on');
clim(ax, ampLimits);
cb = colorbar(ax);
cb.Label.String = 'Amplitude of P_{tot}';
xlabel(ax, 'K_d (\times 10^{-4})');
ylabel(ax, 'A_T (a.u.)');
xlim(ax, xLimits);
ylim(ax, yLimits);

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%%
function component = join_curve_branches(leftBranch, seed, rightBranch)
leftKd = flipud(branch_kd(leftBranch));
leftAT = flipud(branch_at(leftBranch));
leftAmp = flipud(branch_obs_amp(leftBranch));

rightKd = branch_kd(rightBranch);
rightAT = branch_at(rightBranch);
rightAmp = branch_obs_amp(rightBranch);

component = struct();
component.Kd = [leftKd; seed.Parameters(1); rightKd];
component.AT = [leftAT; seed.Parameters(2); rightAT];
component.obsAmp = [leftAmp; seed.obsAmp(1); rightAmp];
end

function component = stitch_nearest_components(lower, upper)
candidates = cell(2, 2);
distances = inf(2, 2);

for i = 1:2
    for j = 1:2
        first = maybe_reverse_component(lower, i == 2);
        second = maybe_reverse_component(upper, j == 2);
        candidates{i, j} = [first, second];
        distances(i, j) = hypot(first.Kd(end) - second.Kd(1), first.AT(end) - second.AT(1));
    end
end

[~, idx] = min(distances(:));
[rowIdx, colIdx] = ind2sub(size(distances), idx);
selected = candidates{rowIdx, colIdx};

component = struct();
component.Kd = [selected(1).Kd; selected(2).Kd];
component.AT = [selected(1).AT; selected(2).AT];
component.obsAmp = [selected(1).obsAmp; selected(2).obsAmp];
end

function component = maybe_reverse_component(component, shouldReverse)
if ~shouldReverse
    return
end

component.Kd = flipud(component.Kd);
component.AT = flipud(component.AT);
component.obsAmp = flipud(component.obsAmp);
end

function anchor = curve_label_anchor(loaded)
branches = {loaded.upperRightBranch, loaded.upperLeftBranch, loaded.lowerRightBranch, loaded.lowerLeftBranch};
for i = 1:numel(branches)
    if isempty(branches{i})
        continue
    end
    params = vertcat(branches{i}.params);
    if ~isempty(params)
        anchor = params(end, :);
        return
    end
end
anchor = loaded.upperSeed.Parameters;
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

function values = branch_obs_amp(branch)
if isempty(branch)
    values = zeros(0, 1);
    return
end
values = arrayfun(@(entry) entry.derived.obsAmp(1), branch).';
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end

function tag = amplitude_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end
