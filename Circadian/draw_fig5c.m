clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
figureFile = fullfile(scriptDir, 'fig5c.png');
markerDataDir = fullfile(scriptDir, 'data', 'fig5b', 'markers');

%% Files
targetPeriod = 24.0;
targetAmplitudes = [0.01, 0.02, 0.04];

markerFiles = cell(1, numel(targetAmplitudes));
for i = 1:numel(targetAmplitudes)
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('fig5b_marker_T%s_A%s.mat', period_tag(targetPeriod), amplitude_tag(targetAmplitudes(i))));
end

%% Plot settings
nColor = 256;
tColor = linspace(0, 1, nColor).';
warmColormap = (1 - tColor) .* [1.00, 0.55, 0.14] + tColor .* [0.98, 0.96, 0.92];
amplitudeLimits = [min(targetAmplitudes), max(targetAmplitudes)];

%% Load marker time series
series = repmat(struct( ...
    't', [], ...
    'Ptot', [], ...
    'Pc', [], ...
    'Pn', [], ...
    'amplitude', []), 1, numel(markerFiles));

for i = 1:numel(markerFiles)
    loaded = load(markerFiles{i}, 'TS', 'obsAmp', 'period');
    t = loaded.TS{1}(:);
    y = loaded.TS{2};
    obs = y(:, 2) + y(:, 3);
    [shiftedT, order] = shift_cycle_to_max_local(t, obs, loaded.period(1));

    series(i).t = shiftedT;
    series(i).Ptot = obs(order);
    series(i).Pc = y(order, 2);
    series(i).Pn = y(order, 3);
    series(i).amplitude = loaded.obsAmp(1);
end

%% Draw figure
fig = figure('Color', 'w', 'Position', [100, 100, 1180, 360]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(series)
    ax = nexttile(tl, i);
    hold(ax, 'on');

    baseColor = value_to_color(series(i).amplitude, amplitudeLimits, warmColormap);
    plot(ax, series(i).t, series(i).Ptot, 'LineWidth', 3.0, 'Color', baseColor);
    plot(ax, series(i).t, series(i).Pc, ':', 'LineWidth', 2.5, 'Color', blend(baseColor, 0.35));
    plot(ax, series(i).t, series(i).Pn, '-.', 'LineWidth', 2.5, 'Color', blend(baseColor, 0.60));

    grid(ax, 'on');
    xlabel(ax, 'Time (hour)');
    ylabel(ax, 'Concentration (a.u.)');
    title(ax, sprintf('T = 24, A = %.2f', targetAmplitudes(i)));
end

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%%
function [shiftedT, order] = shift_cycle_to_max_local(t, obs, period)
obs = obs(:);
[~, peakIdx] = max(obs);
shiftedT = mod(t(:) - t(peakIdx), period);
[shiftedT, order] = sort(shiftedT);
end

function color = value_to_color(value, valueLimits, cmap)
position = 1 + (size(cmap, 1) - 1) * (value - valueLimits(1)) / (valueLimits(2) - valueLimits(1));
position = min(max(position, 1), size(cmap, 1));
lowerIdx = floor(position);
upperIdx = ceil(position);
weight = position - lowerIdx;
color = (1 - weight) * cmap(lowerIdx, :) + weight * cmap(upperIdx, :);
end

function mixed = blend(color, alpha)
mixed = (1 - alpha) * color + alpha * [1, 1, 1];
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end

function tag = amplitude_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end
