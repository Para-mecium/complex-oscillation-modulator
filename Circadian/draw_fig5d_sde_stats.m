clear
clc

%% Paths
disp('%% Paths');
scriptDir = fileparts(mfilename('fullpath'));
sdeDataDir = fullfile(scriptDir, 'data', 'fig5d', 'sde');
repeatDataFile = fullfile(sdeDataDir, 'fig5d_sde_repeat.mat');
psdFigureFile = fullfile(scriptDir, 'fig5d_sde_stats_psd.png');
distributionFigureFile = fullfile(scriptDir, 'fig5d_sde_stats_distribution.png');

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
periodColormap = [
    6, 16, 111
    43, 96, 166
    72, 157, 212
] / 255;
periodEdgeColormap = [
    6, 16, 111
    43, 96, 166
    72, 157, 212
] / 255;
periodValues = [psdStats.targetPeriod];
periodLimits = [min(periodValues), max(periodValues)];
psdFrequencyScale = 24;
psdXAxis = [0.9, 1.1];
psdYAxis = [];
maxPtotScale = 1;
maxPtotXAxis = [0.11, 0.13];
maxPtotYAxis = [0, 0.3];
maxPtotDistributionMode = 'histogram';
maxPtotBinCount = 50;
maxPtotFaceAlphaByPeriod = [0.8, 0.8, 0.8];
maxPtotEdgeAlphaByPeriod = [1,1,1];
maxPtotFitPointCount = 400;
maxPtotFitLineWidth = 2.4;

%% Draw figure
disp('%% Draw figure');
figPsd = figure('Color', 'w');
ax1 = axes(figPsd);
hold(ax1, 'on');
for i = 1:numel(psdStats)
    fprintf('Plot PSD curve %d/%d: T = %.1f\n', ...
        i, numel(psdStats), psdStats(i).targetPeriod);
    lineColor = value_to_color(psdStats(i).targetPeriod, periodLimits, periodColormap);
    scaledFrequency = psdFrequencyScale * psdStats(i).frequency;
    visibleMask = scaledFrequency >= psdXAxis(1) & scaledFrequency <= psdXAxis(2);
    visibleFrequency = scaledFrequency(visibleMask);
    visibleLower = psdStats(i).lower(visibleMask);
    visibleUpper = psdStats(i).upper(visibleMask);
    visibleMean = psdStats(i).mean(visibleMask);
    fill(ax1, [visibleFrequency; flipud(visibleFrequency)], ...
        [visibleLower; flipud(visibleUpper)], ...
        lineColor, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax1, visibleFrequency, visibleMean, 'LineWidth', 1, ...
        'Color', lineColor, ...
        'DisplayName', sprintf('T = %.1f', psdStats(i).targetPeriod));
    xline(ax1, psdFrequencyScale / psdStats(i).targetPeriod, '--', ...
        'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
grid(ax1, 'on');
box(ax1, 'on');
xlabel(ax1, 'Frequency (1/24 hour^{-1})');
ylabel(ax1, 'PSD of P_{tot} (a.u.^2 hour)');
xlim(ax1, psdXAxis);
if isempty(psdYAxis)
    ylim(ax1, resolve_psd_ylim(psdXAxis, psdFrequencyScale, psdStats));
else
    ylim(ax1, psdYAxis);
end
% legend(ax1, 'Location', 'best');

figDistribution = figure('Color', 'w');
ax2 = axes(figDistribution);
hold(ax2, 'on');
histogramEdges = linspace(maxPtotXAxis(1), maxPtotXAxis(2), maxPtotBinCount + 1);
xFit = linspace(maxPtotXAxis(1), maxPtotXAxis(2), maxPtotFitPointCount).';
for i = 1:numel(distributionStats)
    fprintf('Plot distribution %d/%d: T = %.1f\n', ...
        i, numel(distributionStats), distributionStats(i).targetPeriod);
    lineColor = value_to_color(distributionStats(i).targetPeriod, periodLimits, periodColormap);
    edgeColor = value_to_color(distributionStats(i).targetPeriod, periodLimits, periodEdgeColormap);
    faceAlpha = value_to_color(distributionStats(i).targetPeriod, periodLimits, maxPtotFaceAlphaByPeriod.');
    edgeAlpha = value_to_color(distributionStats(i).targetPeriod, periodLimits, maxPtotEdgeAlphaByPeriod.');
    maxPtot = maxPtotScale * distributionStats(i).maxPtot;
    switch maxPtotDistributionMode
        case 'histogram'
            histogram(ax2, maxPtot, ...
                'BinEdges', histogramEdges, ...
                'Normalization', 'probability', ...
                'FaceColor', lineColor, ...
                'FaceAlpha', faceAlpha, ...
                'EdgeColor', edgeColor, ...
                'EdgeAlpha', edgeAlpha, ...
                'LineWidth', 0.5, ...
                'DisplayName', sprintf('T = %.1f', distributionStats(i).targetPeriod));
        case 'fit'
            density = kernel_density(maxPtot, xFit);
            plot(ax2, xFit, density, 'LineWidth', maxPtotFitLineWidth, ...
                'Color', lineColor, ...
                'DisplayName', sprintf('T = %.1f', distributionStats(i).targetPeriod));
        otherwise
            error('Unknown maxPtotDistributionMode: %s', maxPtotDistributionMode);
    end
end
grid(ax2, 'on');
box(ax2, 'on');
xlabel(ax2, 'Maximum (a.u.)');
if strcmp(maxPtotDistributionMode, 'fit')
    ylabel(ax2, 'Density');
else
    ylabel(ax2, 'Probability');
end
xlim(ax2, maxPtotXAxis);
if ~isempty(maxPtotYAxis)
    ylim(ax2, maxPtotYAxis);
end
legend(ax2, 'Location', 'best');

exportgraphics(figPsd, psdFigureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', psdFigureFile);
exportgraphics(figDistribution, distributionFigureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', distributionFigureFile);

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

function density = kernel_density(samples, evaluationPoints)
samples = samples(:);
samples = samples(isfinite(samples));
if isempty(samples)
    density = zeros(size(evaluationPoints));
    return
end

nSample = numel(samples);
bandwidth = 1.06 * std(samples) * nSample ^ (-1 / 5);
if ~(isfinite(bandwidth) && bandwidth > 0)
    bandwidth = max(range(evaluationPoints), eps) / 100;
end

z = (evaluationPoints(:) - samples.') / bandwidth;
density = mean(exp(-0.5 * z.^2), 2) / (bandwidth * sqrt(2 * pi));
end
