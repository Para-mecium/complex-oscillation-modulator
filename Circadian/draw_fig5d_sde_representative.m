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
plotTimeLimits = [0, 100];
plotYLimits = [0.05, 0.25];
periodColormap = [
    0.02, 0.07, 0.45
    0.06, 0.25, 0.62
    0.31, 0.58, 0.78
];
periodValues = [representative.targetPeriod];
periodLimits = [min(periodValues), max(periodValues)];

%% Draw figure
disp('%% Draw figure');
fig = figure();
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(representative)
    fprintf('Draw representative trace %d/%d: T = %.1f\n', ...
        i, numel(representative), representative(i).targetPeriod);
    Ptot = representative(i).X(:, 2) + representative(i).X(:, 3);
    timeMask = representative(i).t >= plotTimeLimits(1) & representative(i).t <= plotTimeLimits(2);
    lineColor = value_to_color(representative(i).targetPeriod, periodLimits, periodColormap);
    plot(ax, representative(i).t(timeMask), Ptot(timeMask), 'LineWidth', 2.2, ...
        'Color', lineColor, ...
        'DisplayName', sprintf('T = %.1f', representative(i).targetPeriod));
end

grid on
box on
xlim(ax, plotTimeLimits);
ylim(ax, plotYLimits);
xticks(ax, [0, 50, 100]);
yticks(ax, [0.05, 0.15, 0.25]);
xlabel(ax, 'Time (hour)');
ylabel(ax, 'Concentration (a.u.)');
legend(ax, 'Location', 'north', 'Box', 'on');

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
