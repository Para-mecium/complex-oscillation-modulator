function result = reproduce_fig3b(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);
cacheFile = fullfile(cfg.io.cacheDir, 'fig3b_data.mat');
data = flexmod.cache_get_or_create(cacheFile, @() build_fig3b_data(cfg));

fig = plot_data(data, cfg);
flexmod.save_figure(fig, fullfile(cfg.io.figureDir, 'fig3b'));

result = struct('data', data, 'figure', fig);
end

function fig = plot_data(data, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 1500, 420]);
tiled = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled, 1, [1, 1]);
hold(ax1, 'on');
scatter(ax1, data.periodGrid(:), data.amplitudeGrid(:), 10, [0.82, 0.82, 0.82], 'filled', ...
    'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.2);
scatter(ax1, data.varyET.period, data.varyET.amplitude, 38, data.varyET.params(:, 2), 'filled');
scatter(ax1, data.varyI1.period, data.varyI1.amplitude, 38, data.varyI1.params(:, 1), 'filled', ...
    'MarkerFaceAlpha', 0.85);
quiver(ax1, data.varyET.period(2), data.varyET.amplitude(2), ...
    data.varyET.period(end - 1) - data.varyET.period(2), ...
    data.varyET.amplitude(end - 1) - data.varyET.amplitude(2), 0, ...
    'Color', [0.8, 0.25, 0.18], 'LineWidth', 1.4, 'MaxHeadSize', 0.4);
quiver(ax1, data.varyI1.period(2), data.varyI1.amplitude(2), ...
    data.varyI1.period(end - 1) - data.varyI1.period(2), ...
    data.varyI1.amplitude(end - 1) - data.varyI1.amplitude(2), 0, ...
    'Color', [0.12, 0.65, 0.7], 'LineWidth', 1.4, 'MaxHeadSize', 0.4);
grid(ax1, 'on');
xlabel(ax1, 'Period (min)');
ylabel(ax1, 'Amplitude (a.u.)');
title(ax1, 'Grid Scan In (T, A) Space');

for i = 1:3
    ax = nexttile(tiled, i + 1);
    orbit = data.seriesOrbits{i};
    if orbit.success
        plot(ax, orbit.t, orbit.y(:, 2), 'LineWidth', 3, 'Color', [0.15 + 0.25 * i, 0.55 + 0.08 * i, 0.85]);
        xLimit = max(200, orbit.period * 2.5);
        yLimit = 1.1 * max(6, orbit.yMax);
    else
        xLimit = 200;
        yLimit = 6;
    end
    grid(ax, 'on');
    xlim(ax, [0, xLimit]);
    ylim(ax, [0, yLimit]);
    xlabel(ax, 'Time (min)');
    ylabel(ax, 'Y (a.u.)');
    title(ax, sprintf('I_1 = %.1f, E_T = %.1f', cfg.grid.fixedI1, data.seriesET(i)));
end
end
