function result = reproduce_fig3c(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);
curveData = build_fig3c_data(cfg);
markData = build_fig3c_marks(cfg);

fig = plot_data(curveData, markData, cfg);
flexmod.save_figure(fig, fullfile(cfg.io.figureDir, 'fig3c'));

data = struct();
data.curves = curveData.curves;
data.curveFiles = curveData.curveFiles;
data.markResults = markData.markResults;
data.markFiles = markData.markFiles;
data.markAmplitudes = markData.markAmplitudes;
result = struct('data', data, 'figure', fig);
end

function fig = plot_data(curveData, markData, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 480]);
tiled = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
bifurcation = base_bifurcation_curve(cfg);

ax1 = nexttile(tiled, 1);
hold(ax1, 'on');
colormap(ax1, cfg.fig3c.warmColormap);
ampLimits = [min(cellfun(@(c) min(c.amplitude), curveData.curves)), ...
    max(cellfun(@(c) max(c.amplitude), curveData.curves))];
ampLimits = normalize_limits(ampLimits);
for i = 1:numel(curveData.curves)
    curve = curveData.curves{i};
    if ~isempty(curve.I1)
        flexmod.plot_gradient_curve(ax1, curve.I1, curve.ET, curve.amplitude, 3);
    end
    text(ax1, curve.I1(end) + 0.02, curve.ET(end), sprintf('%d', curve.targetPeriod), ...
        'FontSize', 11, 'Color', [0.15, 0.35, 0.75]);
end

function limits = normalize_limits(limits)
if ~all(isfinite(limits))
    limits = [0, 1];
elseif limits(1) >= limits(2)
    delta = max(abs(limits(1)) * 0.05, 1e-6);
    limits = [limits(1) - delta, limits(2) + delta];
end
end
for i = 1:numel(markData.markResults)
    mark = markData.markResults{i};
    scatter(ax1, mark.params(1), mark.params(2), 160, mark.measures.amplitude, 'filled', ...
        'MarkerEdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 1.2);
end
for i = 1:numel(bifurcation.visibleHopfCurves)
    curve = bifurcation.visibleHopfCurves{i};
    plot_bifurcation_segments(ax1, curve, cfg.bifurcation.lineColor, ...
        cfg.bifurcation.lineStyle, cfg.bifurcation.lineWidth);
end
for i = 1:numel(bifurcation.visibleLpCurves)
    curve = bifurcation.visibleLpCurves{i};
    plot_bifurcation_segments(ax1, curve, cfg.bifurcation.lineColor, ...
        cfg.bifurcation.upperLineStyle, cfg.bifurcation.lineWidth);
end
grid(ax1, 'on');
clim(ax1, ampLimits);
cb = colorbar(ax1);
cb.Label.String = 'Protein amplitude';
xlabel(ax1, 'I_1 (a.u.)');
ylabel(ax1, 'E_T (a.u.)');
title(ax1, sprintf('Iso-period curves (T = %d highlighted)', cfg.fig3c.markPeriod));
xlim(ax1, cfg.fig3c.I1Range);
ylim(ax1, bifurcation.yRange);

ax2 = nexttile(tiled, 2);
hold(ax2, 'on');
colormap(ax2, cfg.fig3c.warmColormap);
for i = 1:numel(markData.markResults)
    mark = markData.markResults{i};
    color = value_to_color(mark.measures.amplitude, ampLimits, cfg.fig3c.warmColormap);
    plot(ax2, mark.orbit.t, mark.orbit.y(:, 2), 'LineWidth', 3, 'Color', color, ...
        'DisplayName', sprintf('A = %.1f', markData.markAmplitudes(i)));
end
grid(ax2, 'on');
xlabel(ax2, 'Time (min)');
ylabel(ax2, 'Y (a.u.)');
title(ax2, sprintf('Period = %d min', cfg.fig3c.markPeriod));
legend(ax2, 'Location', 'northwest');
end

function color = value_to_color(value, limits, cmap)
alpha = (value - limits(1)) / max(limits(2) - limits(1), eps);
alpha = min(max(alpha, 0), 1);
idx = 1 + round(alpha * (size(cmap, 1) - 1));
color = cmap(idx, :);
end

function plot_bifurcation_segments(ax, curve, color, lineStyle, lineWidth)
if numel(curve.I1) < 2
    return
end

segments = split_curve_segments(curve.I1(:), curve.ET(:));
for i = 1:numel(segments)
    idx = segments{i};
    if numel(idx) < 2
        continue
    end
    plot(ax, curve.I1(idx), curve.ET(idx), ...
        'Color', color, ...
        'LineStyle', lineStyle, ...
        'LineWidth', lineWidth, ...
        'HandleVisibility', 'off');
end
end

function segments = split_curve_segments(x, y)
segments = {};
if numel(x) < 2
    segments{1} = 1:numel(x);
    return
end

dx = abs(diff(x));
dy = abs(diff(y));
dxRef = max(median(dx(dx > 0)), eps);
dyRef = max(median(dy(dy > 0)), eps);
breakMask = dx > 3 * dxRef | dy > max(0.12, 6 * dyRef);

breaks = [0; find(breakMask); numel(x)];
for i = 1:(numel(breaks) - 1)
    segments{end + 1} = (breaks(i) + 1):breaks(i + 1); %#ok<AGROW>
end
end
