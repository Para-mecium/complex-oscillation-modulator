function fig = draw_figS13a()
%% File paths
scriptDir = fileparts(mfilename('fullpath'));
pathFile = fullfile(scriptDir, 'beta_modulation_path.mat');
initFile = fullfile(scriptDir, 'initData.mat');
learnedFile = fullfile(scriptDir, 'beta_target_data.mat');

%% Load data
pathData = load(pathFile);
initData = load(initFile);
learnedData = load(learnedFile);

%% Figure constants
beforeColor = [0.2, 0.45, 0.85];
afterColor = [0.92, 0.74, 0.16];
pathColor = [0.2, 0.2, 0.2];
zoneColor = [0.75, 0.75, 0.75];

figurePosition = [120, 120, 620, 235];
paramBounds = [0, 4, 0, 25];
propertyBounds = [50, 350, 50, 350];
phaseBounds = [150, 350, 50, 350];

paramTicksX = 0:1:4;
paramTicksY = 0:5:25;
propertyTicksX = 50:50:350;
propertyTicksY = 50:50:350;
phaseTicksX = 150:50:350;
phaseTicksY = 50:50:350;

agingZoneS = 200;
agingZoneH = 100;
controlledIdx = [5 6];
lineWidth = 1.6;

xLabelParam = '$\beta_S$';
yLabelParam = '$\beta_H$';
xLabelProperty = 'min S';
yLabelProperty = 'min H';
xLabelPhase = 'S';
yLabelPhase = 'H';
legendLabels = {'Before', 'After'};

%% Prepare plotting arrays
baselineParams = initData.Parameters(controlledIdx);
targetParams = learnedData.Parameters(controlledIdx);
baselineMinSH = initData.varMin(3:4);
targetMinSH = learnedData.varMin(3:4);
paramMatrix = [pathData.params_start(:).'; pathData.curve_params; pathData.params_end(:).'];
paramPath = paramMatrix(:, controlledIdx);
nPathPoint = size(paramPath, 1);
propertyMinS = linspace(baselineMinSH(1), targetMinSH(1), nPathPoint);
propertyMinH = linspace(baselineMinSH(2), targetMinSH(2), nPathPoint);
propertyPath = [propertyMinS(:), propertyMinH(:)];
baselineTSVar = initData.TS{2};
targetTSVar = learnedData.TS{2};

%% Draw figure
fig = figure('Color', 'w', 'Position', figurePosition);
tiled = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled, 1);
plot(paramPath(:, 1), paramPath(:, 2), '-', 'Color', pathColor, 'LineWidth', lineWidth)
hold on
scatter(baselineParams(1), baselineParams(2), 75, beforeColor, 'filled')
scatter(targetParams(1), targetParams(2), 130, afterColor, 'filled', 'p')
box on
grid on
axis(paramBounds)
pbaspect(ax1, [1, 1, 1])
xlim(ax1, paramBounds(1:2))
ylim(ax1, paramBounds(3:4))
xticks(ax1, paramTicksX)
yticks(ax1, paramTicksY)
xlabel(xLabelParam, 'Interpreter', 'latex', 'FontName', 'Arial')
ylabel(yLabelParam, 'Interpreter', 'latex', 'FontName', 'Arial')
set(ax1, 'FontSize', 11, 'Layer', 'top')

ax2 = nexttile(tiled, 2);
patch([propertyBounds(1), agingZoneS, agingZoneS, propertyBounds(1)], ...
    [propertyBounds(3), propertyBounds(3), propertyBounds(4), propertyBounds(4)], zoneColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on
patch([propertyBounds(1), propertyBounds(2), propertyBounds(2), propertyBounds(1)], ...
    [propertyBounds(3), propertyBounds(3), agingZoneH, agingZoneH], zoneColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(propertyPath(:, 1), propertyPath(:, 2), '-', 'Color', pathColor, 'LineWidth', lineWidth)
scatter(baselineMinSH(1), baselineMinSH(2), 75, beforeColor, 'filled')
scatter(targetMinSH(1), targetMinSH(2), 130, afterColor, 'filled', 'p')
box on
grid on
axis(propertyBounds)
pbaspect(ax2, [1, 1, 1])
xlim(ax2, propertyBounds(1:2))
ylim(ax2, propertyBounds(3:4))
xticks(ax2, propertyTicksX)
yticks(ax2, propertyTicksY)
xlabel(xLabelProperty, 'FontName', 'Arial')
ylabel(yLabelProperty, 'FontName', 'Arial')
set(ax2, 'FontSize', 11, 'Layer', 'top')

ax3 = nexttile(tiled, 3);
patch([phaseBounds(1), agingZoneS, agingZoneS, phaseBounds(1)], ...
    [phaseBounds(3), phaseBounds(3), phaseBounds(4), phaseBounds(4)], zoneColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on
patch([phaseBounds(1), phaseBounds(2), phaseBounds(2), phaseBounds(1)], ...
    [phaseBounds(3), phaseBounds(3), agingZoneH, agingZoneH], zoneColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
baselineHandle = plot(baselineTSVar(:, 3), baselineTSVar(:, 4), 'Color', beforeColor, 'LineWidth', lineWidth);
targetHandle = plot(targetTSVar(:, 3), targetTSVar(:, 4), 'Color', afterColor, 'LineWidth', lineWidth);
box on
grid on
axis(phaseBounds)
pbaspect(ax3, [1, 1, 1])
xlim(ax3, phaseBounds(1:2))
ylim(ax3, phaseBounds(3:4))
xticks(ax3, phaseTicksX)
yticks(ax3, phaseTicksY)
xlabel(xLabelPhase, 'FontName', 'Arial')
ylabel(yLabelPhase, 'FontName', 'Arial')
legend([baselineHandle, targetHandle], legendLabels, 'Location', 'best', 'Box', 'off')
set(ax3, 'FontSize', 11, 'Layer', 'top')
end
