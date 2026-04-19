function fig = draw_etd_fig2e()
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'sde_distribution.mat');

data = load(dataFile);
initColor = [0.2, 0.45, 0.85];
modulatedColor = [0.92, 0.74, 0.16];
zoneColor = [0.75, 0.75, 0.75];
agingZoneS = 200;
agingZoneH = 100;

fig = figure('Color', 'w', 'Position', [120, 120, 650, 360]);
tiled = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tiled, 1);
draw_trace_panel(ax, data.init.trajectory.t, data.init.trajectory.S, ...
    initColor, zoneColor, agingZoneS, [0, 400], 'SIR2 protein (a.u.)', 'Init');

ax = nexttile(tiled, 2);
draw_trace_panel(ax, data.modulated.trajectory.t, data.modulated.trajectory.S, ...
    modulatedColor, zoneColor, agingZoneS, [0, 400], 'SIR2 protein (a.u.)', 'Modulated');

ax = nexttile(tiled, 3);
draw_trace_panel(ax, data.init.trajectory.t, data.init.trajectory.H, ...
    initColor, zoneColor, agingZoneH, [0, 600], 'HAP4 protein (a.u.)', '');

ax = nexttile(tiled, 4);
draw_trace_panel(ax, data.modulated.trajectory.t, data.modulated.trajectory.H, ...
    modulatedColor, zoneColor, agingZoneH, [0, 600], 'HAP4 protein (a.u.)', '');
end

function draw_trace_panel(ax, t, values, lineColor, zoneColor, zoneUpper, yLimits, yLabelText, titleText)
meanValues = mean(values, 2);
intervalValues = prctile(values, [10 90], 2);

axes(ax);
patch([t(1), t(end), t(end), t(1)], [yLimits(1), yLimits(1), zoneUpper, zoneUpper], zoneColor, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
patch([t(:); flipud(t(:))], [intervalValues(:, 1); flipud(intervalValues(:, 2))], lineColor, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none');
plot(t, meanValues, 'Color', lineColor, 'LineWidth', 1.5)
box on
grid on
xlim([t(1), t(end)])
ylim(yLimits)
xlabel('Time (hour)', 'FontName', 'Arial')
ylabel(yLabelText, 'FontName', 'Arial')
title(titleText, 'FontName', 'Arial')
set(ax, 'FontSize', 10, 'Layer', 'top')
end
