clear
clc

%% Paths
disp('%% Paths');
scriptDir = fileparts(mfilename('fullpath'));
sdeDataDir = fullfile(scriptDir, 'data', 'fig5d', 'sde');
representativeDataFile = fullfile(sdeDataDir, 'fig5d_sde_representative.mat');
figureFile = fullfile(scriptDir, 'fig5d_sde_representative.png');

%% Load data
disp('%% Load data');
loaded = load(representativeDataFile, ...
    'representative', ...
    'representativeOptions');
representative = loaded.representative;
representativeOptions = loaded.representativeOptions;

%% Plot settings
disp('%% Plot settings');
nColor = 256;
tColor = linspace(0, 1, nColor).';
coolColormap = (1 - tColor) .* [0.05, 0.25, 0.60] + tColor .* [0.92, 0.95, 0.98];
periodValues = [representative.targetPeriod];
periodLimits = [min(periodValues), max(periodValues)];

%% Draw figure
disp('%% Draw figure');
fig = figure('Color', 'w', 'Position', [100, 100, 680, 420]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(representative)
    fprintf('Draw representative trace %d/%d: T = %.1f\n', ...
        i, numel(representative), representative(i).targetPeriod);
    Ptot = representative(i).X(:, 2) + representative(i).X(:, 3);
    lineColor = value_to_color(representative(i).targetPeriod, periodLimits, coolColormap);
    plot(ax, representative(i).t, Ptot, 'LineWidth', 2.2, ...
        'Color', lineColor, ...
        'DisplayName', sprintf('T = %.1f', representative(i).targetPeriod));
end

grid(ax, 'on');
xlabel(ax, 'Time (hour)');
ylabel(ax, 'P_{tot} (a.u.)');
title(ax, 'Fig. 5d2: Representative noisy time series');
legend(ax, 'Location', 'best');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%%
function color = value_to_color(value, valueLimits, cmap)
position = 1 + (size(cmap, 1) - 1) * (value - valueLimits(1)) / (valueLimits(2) - valueLimits(1));
position = min(max(position, 1), size(cmap, 1));
lowerIdx = floor(position);
upperIdx = ceil(position);
weight = position - lowerIdx;
color = (1 - weight) * cmap(lowerIdx, :) + weight * cmap(upperIdx, :);
end
