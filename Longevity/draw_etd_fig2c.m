function fig = draw_etd_fig2c()
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'sde_representative.mat');

if ~isfile(dataFile)
    generate_extended_fig2c_data
end

data = load(dataFile);
initColor = [0.2, 0.45, 0.85];
modulatedColor = [0.92, 0.74, 0.16];
zoneColor = [0.75, 0.75, 0.75];
agingZoneS = 200;
agingZoneH = 100;

fig = figure('Color', 'w', 'Position', [120, 120, 650, 360]);
tiled = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tiled, 1);
draw_trace_panel(ax, data.init.representative.t, data.init.representative.X(:, 3), ...
    initColor, zoneColor, agingZoneS, [0, 400], 'SIR2 protein (a.u.)', 'Init');

ax = nexttile(tiled, 2);
draw_trace_panel(ax, data.modulated.representative.t, data.modulated.representative.X(:, 3), ...
    modulatedColor, zoneColor, agingZoneS, [0, 400], 'SIR2 protein (a.u.)', 'Modulated');

ax = nexttile(tiled, 3);
draw_trace_panel(ax, data.init.representative.t, data.init.representative.X(:, 4), ...
    initColor, zoneColor, agingZoneH, [0, 600], 'HAP4 protein (a.u.)', '');

ax = nexttile(tiled, 4);
draw_trace_panel(ax, data.modulated.representative.t, data.modulated.representative.X(:, 4), ...
    modulatedColor, zoneColor, agingZoneH, [0, 600], 'HAP4 protein (a.u.)', '');
end

function draw_trace_panel(ax, t, values, lineColor, zoneColor, zoneUpper, yLimits, yLabelText, titleText)
axes(ax);
patch([t(1), t(end), t(end), t(1)], [yLimits(1), yLimits(1), zoneUpper, zoneUpper], zoneColor, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(t, values, 'Color', lineColor, 'LineWidth', 1.5)
box on
grid on
xlim([t(1), t(end)])
ylim(yLimits)
xlabel('Time (hour)', 'FontName', 'Arial')
ylabel(yLabelText, 'FontName', 'Arial')
title(titleText, 'FontName', 'Arial')
set(ax, 'FontSize', 10, 'Layer', 'top')
end
