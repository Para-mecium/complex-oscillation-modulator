function fig = draw_longevity_modulation_panels(cache)
sampleStep = 10;
lineWidth = 1.6;
pathColor = [0.2, 0.2, 0.2];
baselineColor = [0.2, 0.45, 0.85];
targetColor = [0.92, 0.74, 0.16];
zoneColor = [0.75, 0.75, 0.75];

profile = get_figure_profile(cache);
paramIdx = reshape(double(cache.meta.controlledIdx), 1, []);
paramNames = reshape(cache.meta.controlledNames, 1, []);
agingZone = cache.meta.agingZone;

pathIdx = 1:sampleStep:size(cache.path.paramMatrix, 1);
if pathIdx(end) ~= size(cache.path.paramMatrix, 1)
    pathIdx = [pathIdx, size(cache.path.paramMatrix, 1)];
end

paramPath = cache.path.paramMatrix(pathIdx, paramIdx);
propertyPath = [cache.path.minS(pathIdx), cache.path.minH(pathIdx)];

fig = figure('Color', 'w', 'Position', profile.figurePosition);
tiled = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled, 1);
plot(paramPath(:, 1), paramPath(:, 2), '-', 'Color', pathColor, 'LineWidth', lineWidth)
hold on
scatter(cache.baseline.Parameters(paramIdx(1)), cache.baseline.Parameters(paramIdx(2)), ...
    75, baselineColor, 'filled')
scatter(cache.target.Parameters(paramIdx(1)), cache.target.Parameters(paramIdx(2)), ...
    130, targetColor, 'filled', 'p')
box on
grid on
axis(profile.paramBounds)
pbaspect(ax1, [1, 1, 1])
xlim(ax1, profile.paramBounds(1:2))
ylim(ax1, profile.paramBounds(3:4))
if ~isempty(profile.paramTicksX)
    xticks(ax1, profile.paramTicksX)
end
if ~isempty(profile.paramTicksY)
    yticks(ax1, profile.paramTicksY)
end
xlabel(format_parameter_label(paramNames{1}), 'Interpreter', 'latex', 'FontName', 'Arial')
ylabel(format_parameter_label(paramNames{2}), 'Interpreter', 'latex', 'FontName', 'Arial')
set(ax1, 'FontSize', 11, 'Layer', 'top')

ax2 = nexttile(tiled, 2);
draw_aging_zone(ax2, profile.propertyBounds(1:2), profile.propertyBounds(3:4), agingZone, zoneColor);
hold on
plot(propertyPath(:, 1), propertyPath(:, 2), '-', 'Color', pathColor, 'LineWidth', lineWidth)
scatter(cache.baseline.varMin(3), cache.baseline.varMin(4), 75, baselineColor, 'filled')
scatter(cache.target.varMin(3), cache.target.varMin(4), 130, targetColor, 'filled', 'p')
box on
grid on
axis(profile.propertyBounds)
pbaspect(ax2, [1, 1, 1])
xlim(ax2, profile.propertyBounds(1:2))
ylim(ax2, profile.propertyBounds(3:4))
if ~isempty(profile.propertyTicksX)
    xticks(ax2, profile.propertyTicksX)
end
if ~isempty(profile.propertyTicksY)
    yticks(ax2, profile.propertyTicksY)
end
xlabel('min S', 'FontName', 'Arial')
ylabel('min H', 'FontName', 'Arial')
set(ax2, 'FontSize', 11, 'Layer', 'top')

ax3 = nexttile(tiled, 3);
draw_aging_zone(ax3, profile.phaseBounds(1:2), profile.phaseBounds(3:4), agingZone, zoneColor);
hold on
baselineHandle = plot(cache.baseline.TS_var(:, 3), cache.baseline.TS_var(:, 4), ...
    'Color', baselineColor, 'LineWidth', lineWidth);
targetHandle = plot(cache.target.TS_var(:, 3), cache.target.TS_var(:, 4), ...
    'Color', targetColor, 'LineWidth', lineWidth);
box on
grid on
axis(profile.phaseBounds)
pbaspect(ax3, [1, 1, 1])
xlim(ax3, profile.phaseBounds(1:2))
ylim(ax3, profile.phaseBounds(3:4))
if ~isempty(profile.phaseTicksX)
    xticks(ax3, profile.phaseTicksX)
end
if ~isempty(profile.phaseTicksY)
    yticks(ax3, profile.phaseTicksY)
end
xlabel('S', 'FontName', 'Arial')
ylabel('H', 'FontName', 'Arial')
legend([baselineHandle, targetHandle], {'Before', 'After'}, 'Location', 'best', 'Box', 'off')
set(ax3, 'FontSize', 11, 'Layer', 'top')
end

function label = format_parameter_label(name)
switch char(name)
    case 'alphaS'
        label = '$\alpha_S$';
    case 'alphaH'
        label = '$\alpha_H$';
    case 'betaS'
        label = '$\beta_S$';
    case 'betaH'
        label = '$\beta_H$';
    otherwise
        label = char(name);
end
end

function profile = get_figure_profile(cache)
figureId = char(cache.meta.figureId);

profile = struct();
profile.figurePosition = [120, 120, 620, 235];
profile.propertyBounds = [50, 350, 50, 400];
profile.propertyTicksX = 50:50:350;
profile.propertyTicksY = 50:50:400;

switch figureId
    case 'extended_fig2b'
        profile.paramBounds = [0, 35, 0, 1800];
        profile.paramTicksX = 0:5:35;
        profile.paramTicksY = 0:200:1800;
        profile.phaseBounds = [50, 400, 50, 400];
        profile.phaseTicksX = 50:50:400;
        profile.phaseTicksY = 50:50:400;
    case 'figS1a'
        profile.paramBounds = [0, 4, 0, 25];
        profile.paramTicksX = 0:1:4;
        profile.paramTicksY = 0:5:25;
        profile.propertyBounds = [50, 350, 50, 350];
        profile.propertyTicksX = 50:50:350;
        profile.propertyTicksY = 50:50:350;
        profile.phaseBounds = [150, 350, 50, 350];
        profile.phaseTicksX = 150:50:350;
        profile.phaseTicksY = 50:50:350;
    otherwise
        profile.paramBounds = infer_bounds(cache.path.paramMatrix(:, reshape(double(cache.meta.controlledIdx), 1, [])));
        profile.paramTicksX = [];
        profile.paramTicksY = [];
        profile.phaseBounds = infer_bounds([cache.baseline.TS_var(:, 3:4); cache.target.TS_var(:, 3:4)]);
        profile.phaseTicksX = [];
        profile.phaseTicksY = [];
end
end

function bounds = infer_bounds(values)
valueMin = min(values, [], 1);
valueMax = max(values, [], 1);
span = valueMax - valueMin;
scale = max(1, max(abs([valueMin; valueMax]), [], 1));
pad = max(0.08 * span, 1e-3 * scale);
bounds = [valueMin(1) - pad(1), valueMax(1) + pad(1), ...
    valueMin(2) - pad(2), valueMax(2) + pad(2)];
end

function draw_aging_zone(ax, xBounds, yBounds, agingZone, zoneColor)
axes(ax);
patch([xBounds(1), agingZone.S, agingZone.S, xBounds(1)], ...
    [yBounds(1), yBounds(1), yBounds(2), yBounds(2)], zoneColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
patch([xBounds(1), xBounds(2), xBounds(2), xBounds(1)], ...
    [yBounds(1), yBounds(1), agingZone.H, agingZone.H], zoneColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
end
