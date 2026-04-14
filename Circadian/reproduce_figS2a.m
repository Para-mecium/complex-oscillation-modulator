function result = reproduce_figS2a(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
curveData = build_figS2a_data(cfg);
markData = build_figS2a_marks(cfg);

fig = plot_data(curveData, markData, cfg);
circadian.save_figure(fig, fullfile(cfg.io.figureDir, 'figS2a'));
data = struct();
data.curves = curveData.curves;
data.curveFiles = curveData.curveFiles;
data.markResults = markData.markResults;
data.markFiles = markData.markFiles;
data.markAmplitude = markData.markAmplitude;
data.markPeriods = markData.markPeriods;
result = struct('data', data, 'figure', fig);
end

function fig = plot_data(curveData, markData, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 980, 420]);
tiled = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

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
xlim(ax1, cfg.plot.kdScale * cfg.figS2a.KdAxis);
ylim(ax1, cfg.figS2a.ATAxis);
title(ax1, 'Fig. S2a1: Iso-amplitude curves');

ax2 = nexttile(tiled, 2);
hold(ax2, 'on');
for i = 1:numel(markData.markResults)
    mark = markData.markResults{i};
    color = value_to_color(mark.measures.period, periodLimits, cfg.plot.periodColormap);
    orbit = circadian.shift_cycle_to_max(mark.orbit);
    plot(ax2, orbit.t, orbit.obs, 'LineWidth', 3, 'Color', color, ...
        'DisplayName', sprintf('T = %.1f', mark.targetPeriod));
end
grid(ax2, 'on');
xlabel(ax2, 'Time (hour)');
ylabel(ax2, 'P_{tot} (a.u.)');
title(ax2, sprintf('Fig. S2a2: Time series (A = %.3f)', markData.markAmplitude));
legend(ax2, 'Location', 'best');
end

function color = value_to_color(value, limits, cmap)
alpha = (value - limits(1)) / max(limits(2) - limits(1), eps);
alpha = min(max(alpha, 0), 1);
idx = 1 + round(alpha * (size(cmap, 1) - 1));
color = cmap(idx, :);
end
