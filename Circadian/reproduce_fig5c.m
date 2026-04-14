function result = reproduce_fig5c(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
cacheFile = fullfile(cfg.io.cacheDir, 'fig5c_data.mat');
data = circadian.cache_get_or_create(cacheFile, @() build_fig5c_data(cfg));

fig = plot_data(data, cfg);
circadian.save_figure(fig, fullfile(cfg.io.figureDir, 'fig5c'));
result = struct('data', data, 'figure', fig);
end

function fig = plot_data(data, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 1180, 360]);
tiled = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
limits = [min(data.markAmplitudes), max(data.markAmplitudes)];

for i = 1:numel(data.series)
    ax = nexttile(tiled, i);
    hold(ax, 'on');
    baseColor = value_to_color(data.series{i}.amplitude, limits, cfg.plot.amplitudeColormap);
    plot(ax, data.series{i}.t, data.series{i}.Ptot, 'LineWidth', 3, 'Color', baseColor);
    plot(ax, data.series{i}.t, data.series{i}.Pc, ':', 'LineWidth', 2.5, 'Color', blend(baseColor, 0.35));
    plot(ax, data.series{i}.t, data.series{i}.Pn, '-.', 'LineWidth', 2.5, 'Color', blend(baseColor, 0.60));
    grid(ax, 'on');
    xlabel(ax, 'Time (hour)');
    ylabel(ax, 'Concentration (a.u.)');
    title(ax, sprintf('T = 24, A = %.2f', data.markAmplitudes(i)));
end
end

function color = value_to_color(value, limits, cmap)
alpha = (value - limits(1)) / max(limits(2) - limits(1), eps);
alpha = min(max(alpha, 0), 1);
idx = 1 + round(alpha * (size(cmap, 1) - 1));
color = cmap(idx, :);
end

function mixed = blend(color, alpha)
mixed = (1 - alpha) * color + alpha * [1, 1, 1];
end
