function fig = draw_etd_fig2e()
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'sde_distribution.mat');

data = load(dataFile);
initColor = [0, 114, 178] / 255;
modulatedColor = [230, 159, 0] / 255;
zoneColor = [0.25, 0.25, 0.25];
agingZoneS = 200;
agingZoneH = 100;
xLimits = [0, 60];
xTicks = [0, 20, 40, 60];
sYLimits = [100, 450];
sYTicks = [100, 200, 300, 400];
hYLimits = [0, 700];
hYTicks = [0, 200, 400, 600];

fig = gobjects(1, 2);

fig(1) = figure('Color', 'w');
tiledInit = tiledlayout(fig(1), 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tiledInit, 1);
draw_trace_panel(ax, data.init.trajectory.t, data.init.trajectory.S, ...
    initColor, zoneColor, agingZoneS, xLimits, xTicks, sYLimits, sYTicks, {'SIR2', 'protein (a.u.)'});

ax = nexttile(tiledInit, 2);
draw_trace_panel(ax, data.init.trajectory.t, data.init.trajectory.H, ...
    initColor, zoneColor, agingZoneH, xLimits, xTicks, hYLimits, hYTicks, {'HAP4', 'protein (a.u.)'});

fig(2) = figure('Color', 'w');
tiledModulated = tiledlayout(fig(2), 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tiledModulated, 1);
draw_trace_panel(ax, data.modulated.trajectory.t, data.modulated.trajectory.S, ...
    modulatedColor, zoneColor, agingZoneS, xLimits, xTicks, sYLimits, sYTicks, {'SIR2', 'protein (a.u.)'});

ax = nexttile(tiledModulated, 2);
draw_trace_panel(ax, data.modulated.trajectory.t, data.modulated.trajectory.H, ...
    modulatedColor, zoneColor, agingZoneH, xLimits, xTicks, hYLimits, hYTicks, {'HAP4', 'protein (a.u.)'});
end

function draw_trace_panel(ax, t, values, lineColor, zoneColor, zoneUpper, xLimits, xTicks, yLimits, yTicks, yLabelText)
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
xlim(xLimits)
xticks(ax, xTicks)
ylim(yLimits)
yticks(ax, yTicks)
xlabel('Time (hour)', 'FontName', 'Arial')
ylabel(yLabelText, 'FontName', 'Arial')
set(ax, 'FontSize', 10, 'Layer', 'top')
end
