function fig = draw_ed_fig3c()
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'sde_representative.mat');

if ~isfile(dataFile)
    build_ed_fig3c_data
end

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
draw_trace_panel(ax, data.init.representative.t, data.init.representative.X(:, 3), ...
    initColor, zoneColor, agingZoneS, xLimits, xTicks, sYLimits, sYTicks, {'SIR2', 'protein (a.u.)'});

ax = nexttile(tiledInit, 2);
draw_trace_panel(ax, data.init.representative.t, data.init.representative.X(:, 4), ...
    initColor, zoneColor, agingZoneH, xLimits, xTicks, hYLimits, hYTicks, {'HAP4', 'protein (a.u.)'});

fig(2) = figure('Color', 'w');
tiledModulated = tiledlayout(fig(2), 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tiledModulated, 1);
draw_trace_panel(ax, data.modulated.representative.t, data.modulated.representative.X(:, 3), ...
    modulatedColor, zoneColor, agingZoneS, xLimits, xTicks, sYLimits, sYTicks, {'SIR2', 'protein (a.u.)'});

ax = nexttile(tiledModulated, 2);
draw_trace_panel(ax, data.modulated.representative.t, data.modulated.representative.X(:, 4), ...
    modulatedColor, zoneColor, agingZoneH, xLimits, xTicks, hYLimits, hYTicks, {'HAP4', 'protein (a.u.)'});
end

function draw_trace_panel(ax, t, values, lineColor, zoneColor, zoneUpper, xLimits, xTicks, yLimits, yTicks, yLabelText)
axes(ax);
patch([t(1), t(end), t(end), t(1)], [yLimits(1), yLimits(1), zoneUpper, zoneUpper], zoneColor, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
plot(t, values, 'Color', lineColor, 'LineWidth', 1.5)
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
