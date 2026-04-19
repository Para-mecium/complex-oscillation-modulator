function fig = draw_etd_fig2d()
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'sde_distribution.mat');

if ~isfile(dataFile)
    generate_extended_fig2d_data
end

data = load(dataFile);
initColor = [0.2, 0.45, 0.85];
modulatedColor = [0.92, 0.74, 0.16];

fig = figure('Color', 'w', 'Position', [120, 120, 360, 420]);
tiled = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tiled, 1);
draw_hist_panel(ax, data.init.distribution.minS, data.modulated.distribution.minS, ...
    initColor, modulatedColor, 'min SIR2 protein (a.u.)');

ax = nexttile(tiled, 2);
draw_hist_panel(ax, data.init.distribution.minH, data.modulated.distribution.minH, ...
    initColor, modulatedColor, 'min HAP4 protein (a.u.)');
legend(ax, {'Init', 'Modulated'}, 'Location', 'best', 'Box', 'off')
end

function draw_hist_panel(ax, beforeValues, afterValues, beforeColor, afterColor, xLabelText)
allValues = [beforeValues(:); afterValues(:)];
edges = linspace(min(allValues), max(allValues), 26);

axes(ax);
histogram(beforeValues, 'BinEdges', edges, 'Normalization', 'probability', ...
    'FaceColor', beforeColor, 'FaceAlpha', 0.55, 'EdgeAlpha', 0.5);
hold on
histogram(afterValues, 'BinEdges', edges, 'Normalization', 'probability', ...
    'FaceColor', afterColor, 'FaceAlpha', 0.55, 'EdgeAlpha', 0.5);
box on
grid on
xlabel(xLabelText, 'FontName', 'Arial')
ylabel('Probability', 'FontName', 'Arial')
set(ax, 'FontSize', 10, 'Layer', 'top')
end
