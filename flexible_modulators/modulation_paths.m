clear
clc

%% Settings
nCurvePoints = 61;
nOrbitMarkers = 10;
orbitPlotTimeWindow = [0, 200];
orbitPlotYLim = [0, 15];

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
dataDir = fullfile(scriptDir, 'data', 'modulation_paths');
saveFile = fullfile(dataDir, 'modulation_paths_data.mat');
figureFile = fullfile(dataDir, 'modulation_paths.png');
orbitFigureFile = fullfile(dataDir, 'modulation_path_orbits.png');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(dataDir);

%% Temperature model and modulation targets
initialParams = [1.8746, 296.8523];
initialState = [1; 0];

startProperty = [3, 80];
endProperty = [6, 40];

u = linspace(0, 1, nCurvePoints).';

straightPropertyPath = [ ...
    3 + 3 * u, ...
    80 - 40 * u];
openCurvePropertyPath = [ ...
    3 + 3 * u - 1.1 * sin(pi * u), ...
    80 - 40 * u + 28 * sin(pi * u)];
closedCurvePropertyPath = [ ...
    3 + 1.0 * sin(2 * pi * u) + 0.55 * (1 - cos(2 * pi * u)), ...
    80 - 14 * sin(2 * pi * u) + 6 * (1 - cos(2 * pi * u))];

obs = [];
M = 50;
PV.name = 'var';
PV.idx = 2;

controlledIdx = [1, 2];
errBound = 1e-6;

continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);
pathContinuationOptions = continuationOptions;
pathContinuationOptions.initialLambdaStep = 1;

orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 400], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

%% Build system and derivatives
run(fullfile(scriptDir, 'System_temp.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(initialParams));

rhs1 = sys{1};
rhs2 = sys{2};
odeFunc = @(t, y, parameter) [ ...
    rhs1(reshape(y, 1, []), parameter); ...
    rhs2(reshape(y, 1, []), parameter)];

%% Initial periodic orbit of the temperature model
initialOrbit = flexmod_forward_orbit(initialParams, initialState, struct('systemName', 'temp'));

initialFMAMState = state(obs, initialParams, initialOrbit.orbit.t, initialOrbit.orbit.y, M, PV);
initialFMAMState.updatePeriod();
initialFMAMState.updateVar2();
initialSolverView = fmam_state_ops.solverViewFromState(initialFMAMState);

%% Move to the common starting point (A = 3, T = 80)
fprintf('Building start point A=%.3g, T=%.3g\n', startProperty(1), startProperty(2));
startPoint = run_solverview_segment(initialSolverView, startProperty, sys, obs, derivatives, ...
    controlledIdx, errBound, continuationOptions, false);

%% Straight path, implemented by dense piecewise-linear continuation
fprintf('Building straight path with %d target points\n', nCurvePoints);
straightPoints = run_property_path(startPoint, straightPropertyPath, sys, obs, derivatives, ...
    controlledIdx, errBound, pathContinuationOptions);
straightEnd = straightPoints(end);
straightPath = pack_path('straight', straightPropertyPath, straightPoints);
straightParameterPath = build_parameter_path(straightPoints);
straightPath.parameterPath = straightParameterPath;

%% Irregular open curve, implemented by dense piecewise-linear continuation
fprintf('Building open curve path with %d target points\n', nCurvePoints);
openCurvePoints = run_property_path(startPoint, openCurvePropertyPath, sys, obs, derivatives, ...
    controlledIdx, errBound, pathContinuationOptions);
openCurvePath = pack_path('openCurve', openCurvePropertyPath, openCurvePoints);
openCurveParameterPath = build_parameter_path(openCurvePoints);
openCurvePath.parameterPath = openCurveParameterPath;

%% Irregular closed curve, implemented by dense piecewise-linear continuation
fprintf('Building closed curve path with %d target points\n', nCurvePoints);
closedCurvePoints = run_property_path(startPoint, closedCurvePropertyPath, sys, obs, derivatives, ...
    controlledIdx, errBound, pathContinuationOptions);
closedCurvePath = pack_path('closedCurve', closedCurvePropertyPath, closedCurvePoints);
closedCurveParameterPath = build_parameter_path(closedCurvePoints);
closedCurvePath.parameterPath = closedCurveParameterPath;

%% Extract sampled periodic orbits
orbitSampleIdx = round(linspace(1, nCurvePoints, nOrbitMarkers));
straightOrbits = extract_path_orbits(straightPoints, orbitSampleIdx, odeFunc, orbitOptions);
openCurveOrbits = extract_path_orbits(openCurvePoints, orbitSampleIdx, odeFunc, orbitOptions);
closedCurveOrbits = extract_path_orbits(closedCurvePoints, orbitSampleIdx, odeFunc, orbitOptions);

%% Save data
save(saveFile, ...
    'nCurvePoints', ...
    'nOrbitMarkers', ...
    'orbitSampleIdx', ...
    'startProperty', ...
    'endProperty', ...
    'straightPath', ...
    'openCurvePath', ...
    'closedCurvePath', ...
    'straightOrbits', ...
    'openCurveOrbits', ...
    'closedCurveOrbits');

fprintf('Saved data: %s\n', saveFile);

%% Draw parameter-space and property-space paths
straightColor = [0.15, 0.15, 0.15];
openColor = [0.2, 0.45, 0.75];
closedColor = [0.75, 0.32, 0.16];
axisFontSize = 20;
labelFontSize = 20;
titleFontSize = 20;
legendFontSize = 20;
orbitAxisFontSize = 12;
orbitLabelFontSize = 13;
startMarkerSize = 180;
endMarkerSize = 230;
sampleMarkerSize = 120;
arrowSampleIndex = ceil(nOrbitMarkers / 2);
sampleU = u(orbitSampleIdx);
arrowU = mean(sampleU(arrowSampleIndex:arrowSampleIndex + 1));

fig = figure('Color', 'w', 'Position', [120, 120, 1260, 560]);
tiled = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axParam = nexttile(tiled);
hold(axParam, 'on');
plot(axParam, straightParameterPath(:, 1), straightParameterPath(:, 2), ...
    '-', 'Color', straightColor, 'LineWidth', 2.5, 'DisplayName', 'direct straight line');
plot(axParam, openCurveParameterPath(:, 1), openCurveParameterPath(:, 2), ...
    '-', 'Color', openColor, 'LineWidth', 2.5, 'DisplayName', 'nonmonotone open path');
plot(axParam, closedCurveParameterPath(:, 1), closedCurveParameterPath(:, 2), ...
    '-', 'Color', closedColor, 'LineWidth', 2.5, 'DisplayName', 'closed path');
straightSampleParams = vertcat(straightPoints(orbitSampleIdx).Parameters);
openCurveSampleParams = vertcat(openCurvePoints(orbitSampleIdx).Parameters);
closedCurveSampleParams = vertcat(closedCurvePoints(orbitSampleIdx).Parameters);
scatter(axParam, straightSampleParams(:, 1), straightSampleParams(:, 2), sampleMarkerSize, straightColor, ...
    'o', 'LineWidth', 1.8, 'HandleVisibility', 'off');
scatter(axParam, openCurveSampleParams(:, 1), openCurveSampleParams(:, 2), sampleMarkerSize, openColor, ...
    's', 'LineWidth', 1.8, 'HandleVisibility', 'off');
scatter(axParam, closedCurveSampleParams(:, 1), closedCurveSampleParams(:, 2), sampleMarkerSize, closedColor, ...
    'd', 'LineWidth', 1.8, 'HandleVisibility', 'off');
add_parameter_arrow(axParam, u, straightPoints, arrowU, straightColor);
add_parameter_arrow(axParam, u, openCurvePoints, arrowU, openColor);
add_parameter_arrow(axParam, u, closedCurvePoints, arrowU, closedColor);
scatter(axParam, startPoint.Parameters(1), startPoint.Parameters(2), startMarkerSize, [0.12, 0.5, 0.85], 'filled', ...
    'DisplayName', 'initial');
scatter(axParam, straightEnd.Parameters(1), straightEnd.Parameters(2), endMarkerSize, [0.98, 0.78, 0.18], 'filled', 'p', ...
    'DisplayName', 'target');
grid(axParam, 'on');
box(axParam, 'on');
xlabel(axParam, 'E_T (a.u.)', 'FontSize', labelFontSize);
ylabel(axParam, 'Temperature (K)', 'FontSize', labelFontSize);
title(axParam, 'Parameter space', 'FontSize', titleFontSize);
legend(axParam, 'Location', 'best', 'FontSize', legendFontSize);
set(axParam, 'FontSize', axisFontSize);
axis(axParam, 'tight');

axProp = nexttile(tiled);
hold(axProp, 'on');
plot(axProp, straightPath.actualPropertyPath(:, 1), straightPath.actualPropertyPath(:, 2), ...
    '-', 'Color', straightColor, 'LineWidth', 2.5, 'DisplayName', 'direct straight line');
plot(axProp, openCurvePath.actualPropertyPath(:, 1), openCurvePath.actualPropertyPath(:, 2), ...
    '-', 'Color', openColor, 'LineWidth', 2.5, 'DisplayName', 'nonmonotone open path');
plot(axProp, closedCurvePath.actualPropertyPath(:, 1), closedCurvePath.actualPropertyPath(:, 2), ...
    '-', 'Color', closedColor, 'LineWidth', 2.5, 'DisplayName', 'closed path');
straightSampleProps = vertcat(straightPoints(orbitSampleIdx).actualProperty);
openCurveSampleProps = vertcat(openCurvePoints(orbitSampleIdx).actualProperty);
closedCurveSampleProps = vertcat(closedCurvePoints(orbitSampleIdx).actualProperty);
scatter(axProp, straightSampleProps(:, 1), straightSampleProps(:, 2), sampleMarkerSize, straightColor, ...
    'o', 'LineWidth', 1.8, 'HandleVisibility', 'off');
scatter(axProp, openCurveSampleProps(:, 1), openCurveSampleProps(:, 2), sampleMarkerSize, openColor, ...
    's', 'LineWidth', 1.8, 'HandleVisibility', 'off');
scatter(axProp, closedCurveSampleProps(:, 1), closedCurveSampleProps(:, 2), sampleMarkerSize, closedColor, ...
    'd', 'LineWidth', 1.8, 'HandleVisibility', 'off');
add_property_arrow(axProp, 'straight', arrowU, straightColor);
add_property_arrow(axProp, 'open', arrowU, openColor);
add_property_arrow(axProp, 'closed', arrowU, closedColor);
scatter(axProp, startProperty(1), startProperty(2), startMarkerSize, [0.12, 0.5, 0.85], 'filled', ...
    'DisplayName', 'initial');
scatter(axProp, endProperty(1), endProperty(2), endMarkerSize, [0.98, 0.78, 0.18], 'filled', 'p', ...
    'DisplayName', 'target');
grid(axProp, 'on');
box(axProp, 'on');
xlabel(axProp, 'Amplitude A', 'FontSize', labelFontSize);
ylabel(axProp, 'Period T (min)', 'FontSize', labelFontSize);
title(axProp, 'Property space', 'FontSize', titleFontSize);
legend(axProp, 'Location', 'best', 'FontSize', legendFontSize);
set(axProp, 'FontSize', axisFontSize);
axis(axProp, 'tight');

exportgraphics(fig, figureFile, 'Resolution', 300);
fprintf('Saved figure: %s\n', figureFile);

%% Draw sampled periodic orbits
orbitFig = figure('Color', 'w', 'Position', [80, 80, 340 * nOrbitMarkers, 780]);
orbitTiled = tiledlayout(orbitFig, 3, nOrbitMarkers, 'TileSpacing', 'compact', 'Padding', 'compact');

draw_orbit_row(orbitTiled, straightOrbits, straightColor, 'direct straight line', orbitPlotTimeWindow, orbitPlotYLim, orbitAxisFontSize, orbitLabelFontSize, false);
draw_orbit_row(orbitTiled, openCurveOrbits, openColor, 'nonmonotone open path', orbitPlotTimeWindow, orbitPlotYLim, orbitAxisFontSize, orbitLabelFontSize, false);
draw_orbit_row(orbitTiled, closedCurveOrbits, closedColor, 'closed path', orbitPlotTimeWindow, orbitPlotYLim, orbitAxisFontSize, orbitLabelFontSize, true);

exportgraphics(orbitFig, orbitFigureFile, 'Resolution', 300);
fprintf('Saved orbit figure: %s\n', orbitFigureFile);

function add_parameter_arrow(ax, u, points, arrowU, color)

params = vertcat(points.Parameters);
du = u(2) - u(1);
p0 = interp1(u, params, arrowU - du);
p1 = interp1(u, params, arrowU + du);
c = interp1(u, params, arrowU);
add_triangle_arrow(ax, c, p1 - p0, color);
end

function add_property_arrow(ax, pathName, arrowU, color)

du = 0.01;
p0 = evaluate_property_path(pathName, arrowU - du);
p1 = evaluate_property_path(pathName, arrowU + du);
c = evaluate_property_path(pathName, arrowU);
add_triangle_arrow(ax, c, p1 - p0, color);
end

function point = evaluate_property_path(pathName, u)

if strcmp(pathName, 'straight')
    point = [3 + 3 * u, 80 - 40 * u];
elseif strcmp(pathName, 'open')
    point = [3 + 3 * u - 1.1 * sin(pi * u), 80 - 40 * u + 28 * sin(pi * u)];
else
    point = [3 + 1.0 * sin(2 * pi * u) + 0.55 * (1 - cos(2 * pi * u)), ...
        80 - 14 * sin(2 * pi * u) + 6 * (1 - cos(2 * pi * u))];
end
end

function add_triangle_arrow(ax, c, direction, color)

limitsX = xlim(ax);
limitsY = ylim(ax);
scale = [diff(limitsX), diff(limitsY)];

center = c ./ scale;
v = direction ./ scale;
v = v / norm(v);
w = [-v(2), v(1)];

arrowLength = 0.02;
arrowWidth = 0.013;
triangle = [ ...
    center + arrowLength * v; ...
    center - 0.65 * arrowLength * v + arrowWidth * w; ...
    center - 0.65 * arrowLength * v - arrowWidth * w] .* scale;

patch(ax, triangle(:, 1), triangle(:, 2), color, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end

function orbits = extract_path_orbits(points, sampleIdx, odeFunc, orbitOptions)

orbits = repmat(struct('sampleIndex', [], 'Parameters', [], 'actualProperty', [], ...
    't', [], 'y', [], 'period', []), 1, numel(sampleIdx));
for i = 1:numel(sampleIdx)
    point = points(sampleIdx(i));
    y0 = point.TS{2}(1, :).';
    poResult = extract_periodic_orbit(odeFunc, y0, point.Parameters, orbitOptions);

    orbits(i).sampleIndex = sampleIdx(i);
    orbits(i).Parameters = point.Parameters;
    orbits(i).actualProperty = point.actualProperty;
    orbits(i).t = poResult.orbit_t(:);
    orbits(i).y = poResult.orbit_y;
    orbits(i).period = poResult.period;
end
end

function draw_orbit_row(tiled, orbits, color, rowLabel, timeWindow, yRange, axisFontSize, labelFontSize, showXLabel)

for i = 1:numel(orbits)
    ax = nexttile(tiled);
    hold(ax, 'on');
    [tPlot, yPlot] = repeat_orbit(orbits(i), timeWindow(2));
    plot(ax, tPlot, yPlot, 'Color', color, 'LineWidth', 2.0);
    xlim(ax, timeWindow);
    ylim(ax, yRange);
    grid(ax, 'on');
    box(ax, 'on');
    title(ax, sprintf('sample %d', i), 'FontSize', labelFontSize);
    if i == 1
        ylabel(ax, {rowLabel; 'Y'}, 'FontSize', labelFontSize);
    else
        ax.YTickLabel = [];
    end
    if showXLabel
        xlabel(ax, 'Time (min)', 'FontSize', labelFontSize);
    else
        ax.XTickLabel = [];
    end
    set(ax, 'FontSize', axisFontSize);
end
end

function [tPlot, yPlot] = repeat_orbit(orbit, tMax)

t = orbit.t - orbit.t(1);
y = orbit.y(:, 2);
tPlot = [];
yPlot = [];
for k = 0:ceil(tMax / orbit.period)
    tPlot = [tPlot; t + k * orbit.period]; %#ok<AGROW>
    yPlot = [yPlot; y]; %#ok<AGROW>
end
idx = tPlot <= tMax;
tPlot = tPlot(idx);
yPlot = yPlot(idx);
end

function points = run_property_path(startPoint, propertyPath, sys, obs, derivatives, ...
    controlledIdx, errBound, continuationOptions)

points = repmat(startPoint, 1, size(propertyPath, 1));
for i = 2:size(propertyPath, 1)
    fprintf('  segment %d/%d: A=%.4g, T=%.4g\n', ...
        i - 1, size(propertyPath, 1) - 1, propertyPath(i, 1), propertyPath(i, 2));
    points(i) = run_solverview_segment(points(i - 1).solverView, propertyPath(i, :), sys, obs, derivatives, ...
        controlledIdx, errBound, continuationOptions, true);
end
end

function point = run_solverview_segment(seedSolverView, targetProperty, sys, obs, derivatives, ...
    controlledIdx, errBound, continuationOptions, needLog)

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = targetProperty(2) / (2 * pi);

itemsPer(2).prop = 'varAmp';
itemsPer(2).idx = 2;
itemsPer(2).target = targetProperty(1);

task = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
task.psiUpdateMode = true;
task.refreshPsiModeReferences();
task.needLog = needLog;
task.verbose = false;

task.fit();
task.step();
task.errBound = 1e-12;
task.fit();

point = extract_solverview_point_from_task(task, targetProperty);
end

function point = extract_solverview_point_from_task(task, targetProperty)

solverView = task.exportSolverView();
derivedView = task.exportDerivedView();

point = struct();
point.Parameters = reshape(solverView.params, 1, []);
point.TS = {derivedView.t, derivedView.TS_var};
point.period = derivedView.period;
point.varAmp = reshape(derivedView.varAmp, 1, []);
point.varMax = reshape(derivedView.varMax, 1, []);
point.varMin = reshape(derivedView.varMin, 1, []);
point.targetProperty = targetProperty;
point.actualProperty = [point.varAmp(2), point.period];
point.solverView = solverView;
point.derivedView = derivedView;
point.logs = task.logs;
end

function path = pack_path(name, targetPropertyPath, points)

path = struct();
path.name = name;
path.targetPropertyPath = targetPropertyPath;
path.points = points;
path.actualPropertyPath = vertcat(points.actualProperty);
end

function parameterPath = build_parameter_path(points)

parameterPath = points(1).Parameters;
for i = 2:numel(points)
    segmentPath = vertcat(points(i).logs.params);
    if norm(segmentPath(1, :) - parameterPath(end, :), inf) < 1e-10
        segmentPath(1, :) = [];
    end
    parameterPath = [parameterPath; segmentPath]; %#ok<AGROW>
end
end
