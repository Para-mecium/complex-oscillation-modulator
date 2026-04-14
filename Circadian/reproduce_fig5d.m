function result = reproduce_fig5d(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
sdeSigma = resolve_plot_sigma(cfg.fig5d.sde);
curveData = build_fig5d_data(cfg);
markData = build_fig5d_marks(cfg);
sdeData = build_fig5d_sde_data(cfg, markData);

fig = plot_data(curveData, markData, sdeData, cfg);
figureBase = fullfile(cfg.io.figureDir, ...
    sprintf('fig5d_sigma_%s', circadian.format_sigma_tag(sdeSigma)));
circadian.save_figure(fig, figureBase);
data = struct();
data.curves = curveData.curves;
data.curveFiles = curveData.curveFiles;
data.markResults = markData.markResults;
data.markFiles = markData.markFiles;
data.markPeriods = markData.markPeriods;
data.sde = sdeData;
data.figureBase = figureBase;
result = struct('data', data, 'figure', fig);
end

function fig = plot_data(curveData, markData, sdeData, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 1100, 760]);
tiled = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled, 1);
hold(ax1, 'on');
colormap(ax1, cfg.plot.periodColormap);
periodLimits = [min(cellfun(@(c) min(c.period), curveData.curves)), ...
    max(cellfun(@(c) max(c.period), curveData.curves))];
for i = 1:numel(curveData.curves)
    curve = curveData.curves{i};
    if isempty(curve.Kd)
        continue
    end
    circadian.plot_gradient_curve(ax1, cfg.plot.kdScale * curve.Kd, curve.AT, curve.period, 3);
end
for i = 1:numel(markData.markResults)
    mark = markData.markResults{i};
    scatter(ax1, cfg.plot.kdScale * mark.params(1), mark.params(2), 120, mark.measures.period, ...
        'filled', 'MarkerEdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);
end
grid(ax1, 'on');
clim(ax1, periodLimits);
cb = colorbar(ax1);
cb.Label.String = 'Period (hour)';
xlabel(ax1, 'K_d (\times 10^{-4})');
ylabel(ax1, 'A_T (a.u.)');
xlim(ax1, cfg.plot.kdScale * cfg.fig5d.KdAxis);
ylim(ax1, cfg.fig5d.ATAxis);
title(ax1, 'Fig. 5d1: Iso-maximum curves');

ax2 = nexttile(tiled, 2);
hold(ax2, 'on');
for i = 1:numel(sdeData.representative)
    representative = sdeData.representative(i);
    color = value_to_color(representative.targetPeriod, periodLimits, cfg.plot.periodColormap);
    plot(ax2, representative.t, representative.Ptot, 'LineWidth', 2.2, 'Color', color, ...
        'DisplayName', sprintf('T = %.1f', representative.targetPeriod));
end
grid(ax2, 'on');
xlabel(ax2, 'Time (hour)');
ylabel(ax2, 'P_{tot} (a.u.)');
title(ax2, 'Fig. 5d2: Representative noisy time series');
legend(ax2, 'Location', 'best');

ax3 = nexttile(tiled, 3);
hold(ax3, 'on');
psdCfg = cfg.fig5d.sde;
for i = 1:numel(sdeData.psd)
    stats = sdeData.psd(i);
    color = value_to_color(stats.targetPeriod, periodLimits, cfg.plot.periodColormap);
    freqScaled = psdCfg.psdFrequencyScale * stats.frequency;
    fill(ax3, [freqScaled; flipud(freqScaled)], ...
        [stats.lower; flipud(stats.upper)], color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax3, freqScaled, stats.mean, 'LineWidth', 2.2, 'Color', color, ...
        'DisplayName', sprintf('T = %.1f', stats.targetPeriod));
end
grid(ax3, 'on');
xlabel(ax3, 'Frequency (1/24 hour^{-1})');
ylabel(ax3, 'PSD');
title(ax3, sprintf('Fig. 5d3: Mean PSD with %d-%d%% band', ...
    sdeData.meta.psdBandPercentiles(1), sdeData.meta.psdBandPercentiles(2)));
legend(ax3, 'Location', 'best');
xlim(ax3, psdCfg.psdFrequencyAxis);
ylim(ax3, resolve_psd_ylim(psdCfg, sdeData.psd));
ax3.YAxis.Exponent = psdCfg.psdYAxisExponent;
for targetPeriod = reshape(markData.markPeriods, 1, [])
    xline(ax3, psdCfg.psdFrequencyScale / targetPeriod, '--', ...
        'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

ax4 = nexttile(tiled, 4);
hold(ax4, 'on');
allValues = vertcat(sdeData.distribution.maxPtot);
valueMin = min(allValues);
valueMax = max(allValues);
if abs(valueMax - valueMin) < eps(max(abs([valueMin, valueMax]), 1))
    valueMax = valueMin + 1e-6;
end
edges = linspace(valueMin, valueMax, 26);
for i = 1:numel(sdeData.distribution)
    stats = sdeData.distribution(i);
    color = value_to_color(stats.targetPeriod, periodLimits, cfg.plot.periodColormap);
    histogram(ax4, stats.maxPtot, 'BinEdges', edges, 'Normalization', 'probability', ...
        'FaceColor', color, 'FaceAlpha', 0.35, 'EdgeAlpha', 0.45, ...
        'DisplayName', sprintf('T = %.1f', stats.targetPeriod));
end
grid(ax4, 'on');
xlabel(ax4, 'max P_{tot} (a.u.)');
ylabel(ax4, 'Probability');
title(ax4, sprintf('Fig. 5d4: Maximum distribution (n = %d)', sdeData.meta.nReplicates));
legend(ax4, 'Location', 'best');
end

function color = value_to_color(value, limits, cmap)
alpha = (value - limits(1)) / max(limits(2) - limits(1), eps);
alpha = min(max(alpha, 0), 1);
idx = 1 + round(alpha * (size(cmap, 1) - 1));
color = cmap(idx, :);
end

function limits = resolve_psd_ylim(psdCfg, psdStats)
if ~isempty(psdCfg.psdYAxis)
    limits = psdCfg.psdYAxis;
    return
end

visibleMax = 0;
for i = 1:numel(psdStats)
    freqScaled = psdCfg.psdFrequencyScale * psdStats(i).frequency;
    visibleMask = freqScaled >= psdCfg.psdFrequencyAxis(1) & freqScaled <= psdCfg.psdFrequencyAxis(2);
    if ~any(visibleMask)
        continue
    end
    visibleMax = max(visibleMax, max(psdStats(i).upper(visibleMask)));
end

if ~(isfinite(visibleMax) && visibleMax > 0)
    visibleMax = 1;
end

limits = [0, 1.1 * visibleMax];
end

function sigma = resolve_plot_sigma(sdeCfg)
if isfield(sdeCfg, 'plotSigma') && ~isempty(sdeCfg.plotSigma)
    sigma = double(sdeCfg.plotSigma);
else
    sigma = double(sdeCfg.sigma);
end
end
