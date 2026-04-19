clear
clc

%% Paths
disp('%% Paths');
scriptDir = fileparts(mfilename('fullpath'));
sdeDataDir = fullfile(scriptDir, 'data', 'fig5d', 'sde');
repeatDataFile = fullfile(sdeDataDir, 'fig5d_sde_repeat.mat');
figureFile = fullfile(scriptDir, 'fig5d_sde_stats.png');

%% Load data
disp('%% Load data');
loaded = load(repeatDataFile, ...
    'distributionStats', ...
    'psdBandPercentiles', ...
    'psdStats', ...
    'repeatSeeds');
psdStats = loaded.psdStats;
distributionStats = loaded.distributionStats;
psdBandPercentiles = loaded.psdBandPercentiles;
repeatSeeds = loaded.repeatSeeds;

%% Plot settings
disp('%% Plot settings');
nColor = 256;
tColor = linspace(0, 1, nColor).';
coolColormap = (1 - tColor) .* [0.05, 0.25, 0.60] + tColor .* [0.92, 0.95, 0.98];
periodValues = [psdStats.targetPeriod];
periodLimits = [min(periodValues), max(periodValues)];
psdFrequencyScale = 24;
psdFrequencyAxis = [0.9, 1.1];

%% Draw figure
disp('%% Draw figure');
fig = figure('Color', 'w', 'Position', [100, 100, 760, 420]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on');
for i = 1:numel(psdStats)
    fprintf('Plot PSD curve %d/%d: T = %.1f\n', ...
        i, numel(psdStats), psdStats(i).targetPeriod);
    lineColor = value_to_color(psdStats(i).targetPeriod, periodLimits, coolColormap);
    scaledFrequency = psdFrequencyScale * psdStats(i).frequency;
    visibleMask = scaledFrequency >= psdFrequencyAxis(1) & scaledFrequency <= psdFrequencyAxis(2);
    visibleFrequency = scaledFrequency(visibleMask);
    visibleLower = psdStats(i).lower(visibleMask);
    visibleUpper = psdStats(i).upper(visibleMask);
    visibleMean = psdStats(i).mean(visibleMask);
    fill(ax1, [visibleFrequency; flipud(visibleFrequency)], ...
        [visibleLower; flipud(visibleUpper)], ...
        lineColor, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax1, visibleFrequency, visibleMean, 'LineWidth', 2.2, ...
        'Color', lineColor, ...
        'DisplayName', sprintf('T = %.1f', psdStats(i).targetPeriod));
    xline(ax1, psdFrequencyScale / psdStats(i).targetPeriod, '--', ...
        'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
grid(ax1, 'on');
xlabel(ax1, 'Frequency (1/24 hour^{-1})');
ylabel(ax1, 'PSD');
xlim(ax1, psdFrequencyAxis);
ylim(ax1, resolve_psd_ylim(psdFrequencyAxis, psdFrequencyScale, psdStats));
title(ax1, sprintf('Fig. 5d3: Mean PSD with %d-%d%% band', ...
    psdBandPercentiles(1), psdBandPercentiles(2)));
legend(ax1, 'Location', 'best');

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
allMaxima = vertcat(distributionStats.maxPtot);
histogramEdges = linspace(min(allMaxima), max(allMaxima), 26);
for i = 1:numel(distributionStats)
    fprintf('Plot distribution %d/%d: T = %.1f\n', ...
        i, numel(distributionStats), distributionStats(i).targetPeriod);
    lineColor = value_to_color(distributionStats(i).targetPeriod, periodLimits, coolColormap);
    histogram(ax2, distributionStats(i).maxPtot, ...
        'BinEdges', histogramEdges, ...
        'Normalization', 'probability', ...
        'FaceColor', lineColor, ...
        'FaceAlpha', 0.35, ...
        'EdgeAlpha', 0.45, ...
        'DisplayName', sprintf('T = %.1f', distributionStats(i).targetPeriod));
end
grid(ax2, 'on');
xlabel(ax2, 'max P_{tot} (a.u.)');
ylabel(ax2, 'Probability');
title(ax2, sprintf('Fig. 5d4: Maximum distribution (n = %d)', numel(repeatSeeds)));
legend(ax2, 'Location', 'best');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

function limits = resolve_psd_ylim(psdFrequencyAxis, psdFrequencyScale, psdStats)
visibleMax = 0;

for i = 1:numel(psdStats)
    scaledFrequency = psdFrequencyScale * psdStats(i).frequency;
    visibleMask = scaledFrequency >= psdFrequencyAxis(1) & scaledFrequency <= psdFrequencyAxis(2);
    if any(visibleMask)
        visibleMax = max(visibleMax, max(psdStats(i).upper(visibleMask)));
    end
end

if ~(isfinite(visibleMax) && visibleMax > 0)
    visibleMax = 1;
end

limits = [0, 1.1 * visibleMax];
end

function color = value_to_color(value, valueLimits, cmap)
position = 1 + (size(cmap, 1) - 1) * (value - valueLimits(1)) / (valueLimits(2) - valueLimits(1));
position = min(max(position, 1), size(cmap, 1));
lowerIdx = floor(position);
upperIdx = ceil(position);
weight = position - lowerIdx;
color = (1 - weight) * cmap(lowerIdx, :) + weight * cmap(upperIdx, :);
end
