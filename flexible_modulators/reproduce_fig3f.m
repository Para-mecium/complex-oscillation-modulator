function result = reproduce_fig3f(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);
cacheFile = fullfile(cfg.io.cacheDir, 'fig3f_data.mat');
data = flexmod.cache_get_or_create(cacheFile, @() build_fig3f_data(cfg));

fig = plot_data(data);
flexmod.save_figure(fig, fullfile(cfg.io.figureDir, 'fig3f'));

result = struct('data', data, 'figure', fig);
end

function fig = plot_data(data)
fig = figure('Color', 'w', 'Position', [100, 80, 1380, 760]);
tiled = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

axPath = nexttile(tiled, [3, 1]);
hold(axPath, 'on');
plot(axPath, data.directPath.path.paramMatrix(:, 1), data.directPath.path.paramMatrix(:, 2), ...
    'Color', [0.42, 0.64, 0.1], 'LineWidth', 3, 'DisplayName', 'Direct modulation path');
plot(axPath, data.orthPeriod.path.paramMatrix(:, 1), data.orthPeriod.path.paramMatrix(:, 2), ...
    'Color', [0.98, 0.72, 0.25], 'LineWidth', 3, 'DisplayName', 'Orthogonal period path');
plot(axPath, data.orthAmplitude.path.paramMatrix(:, 1), data.orthAmplitude.path.paramMatrix(:, 2), ...
    'Color', [0.42, 0.7, 0.92], 'LineWidth', 3, 'DisplayName', 'Orthogonal amplitude path');
scatter(axPath, data.startPoint.params(1), data.startPoint.params(2), 220, [0.12, 0.5, 0.85], 'filled');
scatter(axPath, data.directMid.params(1), data.directMid.params(2), 220, [0.85, 0.85, 0.85], 'filled', 'd');
scatter(axPath, data.orthPeriod.path.paramMatrix(end, 1), data.orthPeriod.path.paramMatrix(end, 2), 220, [0.85, 0.85, 0.85], 'filled', 's');
scatter(axPath, data.directPath.params(1), data.directPath.params(2), 320, [0.98, 0.78, 0.18], 'filled', 'p');
grid(axPath, 'on');
xlabel(axPath, 'E_T (a.u.)');
ylabel(axPath, 'Temperature (K)');
legend(axPath, 'Location', 'northwest');
title(axPath, 'Direct vs orthogonal implementations');

plot_orbit_tile(nexttile(tiled), data.startPoint.orbit, [0.12, 0.5, 0.85], 'Direct start');
plot_orbit_tile(nexttile(tiled), data.startPoint.orbit, [0.12, 0.5, 0.85], 'Orthogonal start');
plot_orbit_tile(nexttile(tiled), data.directMid.orbit, [0.86, 0.86, 0.86], 'Direct middle');
plot_orbit_tile(nexttile(tiled), data.orthPeriod.orbit, [0.86, 0.86, 0.86], 'Orthogonal middle');
plot_orbit_tile(nexttile(tiled), data.directPath.orbit, [0.98, 0.78, 0.18], 'Direct end');
plot_orbit_tile(nexttile(tiled), data.orthAmplitude.orbit, [0.98, 0.78, 0.18], 'Orthogonal end');
end

function plot_orbit_tile(ax, orbit, color, titleText)
plot(ax, orbit.t, orbit.y(:, 2), 'LineWidth', 3, 'Color', color);
grid(ax, 'on');
xlabel(ax, 'Time (min)');
ylabel(ax, 'Y (a.u.)');
title(ax, titleText);
end
