function data = build_fig5d_data(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

maxima = cfg.fig5d.maxima;
ATseed = cfg.fig5d.ATseed;
Kdseed = cfg.fig5d.Kdseed;
assert(numel(maxima) == numel(ATseed), ...
    'cfg.fig5d.maxima and cfg.fig5d.ATseed must align.');
assert(numel(maxima) == numel(Kdseed), ...
    'cfg.fig5d.maxima and cfg.fig5d.Kdseed must align.');

curves = cell(1, numel(maxima));
curveFiles = cell(1, numel(maxima));
for i = 1:numel(maxima)
    curveFiles{i} = circadian.curve_cache_file(cfg, 'fig5d', 'maximum', maxima(i));
    curves{i} = circadian.cache_get_or_create(curveFiles{i}, ...
        @() build_iso_maximum_curve(cfg, maxima(i), Kdseed(i)));
end

data = struct();
data.curves = curves;
data.curveFiles = curveFiles;
end

function curve = build_iso_maximum_curve(cfg, targetMaximum, targetseed)
seed = run_modulation_task(struct( ...
    'startParams', cfg.clock.initialParams, ...
    'goalOrder', {{'maximum', 'Kd'}}, ...
    'goals', struct('maximum', targetMaximum, 'Kd', targetseed), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', 1)));

left = run_modulation_task(struct( ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'maximum', 'Kd'}}, ...
    'goals', struct('maximum', targetMaximum, 'Kd', cfg.fig5d.KdRange(1)), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

right = run_modulation_task(struct( ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'maximum', 'Kd'}}, ...
    'goals', struct('maximum', targetMaximum, 'Kd', cfg.fig5d.KdRange(2)), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

leftPath = circadian.trim_path(left.path, 'iso_maximum', cfg);
rightPath = circadian.trim_path(right.path, 'iso_maximum', cfg);
curve = combine_branches(leftPath, rightPath);
curve.targetMaximum = targetMaximum;
curve.figureId = 'fig5d';
curve.featureKind = 'maximum';
curve.targetValue = targetMaximum;
end

function curve = combine_branches(pathA, pathB)
iKd = circadian.parameter_index(pathA, 'Kd');
iAT = circadian.parameter_index(pathA, 'AT');

orderA = sortrows([(1:size(pathA.paramMatrix, 1)).', pathA.paramMatrix(:, iKd)], 2);
orderB = sortrows([(1:size(pathB.paramMatrix, 1)).', pathB.paramMatrix(:, iKd)], 2);
idxA = orderA(:, 1);
idxB = orderB(:, 1);

matrixA = pathA.paramMatrix(idxA, :);
matrixB = pathB.paramMatrix(idxB, :);
periodA = pathA.period(idxA);
periodB = pathB.period(idxB);
maxA = pathA.obsMax(idxA);
maxB = pathB.obsMax(idxB);
ampA = pathA.obsAmplitude(idxA);
ampB = pathB.obsAmplitude(idxB);
condA = pathA.directConditionEstimate(idxA);
condB = pathB.directConditionEstimate(idxB);

if ~isempty(matrixA) && ~isempty(matrixB) && norm(matrixA(end, :) - matrixB(1, :), inf) < 1e-10
    matrixB(1, :) = [];
    periodB(1) = [];
    maxB(1) = [];
    ampB(1) = [];
    condB(1) = [];
end

curve = struct();
curve.paramNames = pathA.paramNames;
curve.paramMatrix = [matrixA; matrixB];
curve.Kd = curve.paramMatrix(:, iKd);
curve.AT = curve.paramMatrix(:, iAT);
curve.period = [periodA; periodB];
curve.obsMax = [maxA; maxB];
curve.obsAmplitude = [ampA; ampB];
curve.directConditionEstimate = [condA; condB];
curve.stopReason = {pathA.stopReason, pathB.stopReason};
curve.segments = { ...
    struct('name', 'left', 'Kd', matrixA(:, iKd), 'AT', matrixA(:, iAT), ...
        'period', periodA, 'obsMax', maxA, 'obsAmplitude', ampA, 'directConditionEstimate', condA), ...
    struct('name', 'right', 'Kd', matrixB(:, iKd), 'AT', matrixB(:, iAT), ...
        'period', periodB, 'obsMax', maxB, 'obsAmplitude', ampB, 'directConditionEstimate', condB)};
end
