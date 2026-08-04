function fig = draw_ed_fig3d()
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'sde_distribution.mat');

if ~isfile(dataFile)
    build_ed_fig3d_data
end

data = load(dataFile);
initColor = [0.2, 0.45, 0.85];
modulatedColor = [0.92, 0.74, 0.16];
zoneColor = [0.25, 0.25, 0.25];
zoneAlpha = 0.25;
sZoneUpper = 200;
hZoneUpper = 100;
sXLimits = [0, 300];
sXTicks = [0, 100, 200, 300];
sBinCount = 75;
sYLimits = [0, 0.25];
sYTicks = 0:0.05:0.25;
hXLimits = [0, 300];
hXTicks = [0, 100, 200, 300];
hBinCount = 75;
hYLimits = [0, 0.25];
hYTicks = 0:0.05:0.25;

fig = gobjects(1, 2);

fig(1) = figure('Color', 'w');
ax = axes(fig(1));
draw_hist_panel(ax, data.init.distribution.minS, data.modulated.distribution.minS, ...
    initColor, modulatedColor, zoneColor, zoneAlpha, sZoneUpper, ...
    sXLimits, sXTicks, sBinCount, sYLimits, sYTicks, 'min SIR2 protein (a.u.)');

fig(2) = figure('Color', 'w');
ax = axes(fig(2));
draw_hist_panel(ax, data.init.distribution.minH, data.modulated.distribution.minH, ...
    initColor, modulatedColor, zoneColor, zoneAlpha, hZoneUpper, ...
    hXLimits, hXTicks, hBinCount, hYLimits, hYTicks, 'min HAP4 protein (a.u.)');
% legend(ax, {'Init', 'Modulated'}, 'Location', 'best', 'Box', 'off')
end

function draw_hist_panel(ax, beforeValues, afterValues, beforeColor, afterColor, zoneColor, zoneAlpha, zoneUpper, xLimits, xTicks, binCount, yLimits, yTicks, xLabelText)
edges = linspace(xLimits(1), xLimits(2), binCount + 1);
axes(ax);
patch([xLimits(1), zoneUpper, zoneUpper, xLimits(1)], ...
    [yLimits(1), yLimits(1), yLimits(2), yLimits(2)], zoneColor, ...
    'FaceAlpha', zoneAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on
histogram(beforeValues, 'BinEdges', edges, 'Normalization', 'probability', ...
    'FaceColor', beforeColor, 'FaceAlpha', 0.55, 'EdgeAlpha', 0.5);
histogram(afterValues, 'BinEdges', edges, 'Normalization', 'probability', ...
    'FaceColor', afterColor, 'FaceAlpha', 0.55, 'EdgeAlpha', 0.5);
box on
grid on
xlim(xLimits)
xticks(ax, xTicks)
ylim(yLimits)
yticks(ax, yTicks)
xlabel(xLabelText, 'FontName', 'Arial')
ylabel('Probability', 'FontName', 'Arial')
set(ax, 'FontSize', 10, 'Layer', 'top')
end
