function data = build_figS2a_data(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

amplitudes = cfg.figS2a.amplitudes;
curves = cell(1, numel(amplitudes));
curveFiles = cell(1, numel(amplitudes));
for i = 1:numel(amplitudes)
    curveFiles{i} = circadian.curve_cache_file(cfg, 'figS2a', 'amplitude', amplitudes(i));
    curves{i} = circadian.cache_get_or_create(curveFiles{i}, ...
        @() build_iso_amplitude_curve(cfg, amplitudes(i)));
end

data = struct();
data.curves = curves;
data.curveFiles = curveFiles;
end

function curve = build_iso_amplitude_curve(cfg, targetAmplitude)
seed = run_circadian_fmam_task(struct( ...
    'startParams', cfg.clock.initialParams, ...
    'goalOrder', {{'amplitude', 'period'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'period', cfg.figS2a.centerPeriod), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', 1)));

left = run_circadian_fmam_task(struct( ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'amplitude', 'Kd'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'Kd', cfg.figS2a.KdRange(1)), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

right = run_circadian_fmam_task(struct( ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'amplitude', 'Kd'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'Kd', cfg.figS2a.KdRange(2)), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

leftPath = circadian.trim_path(left.path, 'iso_amplitude', cfg);
rightPath = circadian.trim_path(right.path, 'iso_amplitude', cfg);
curve = combine_branches(leftPath, rightPath);
curve.targetAmplitude = targetAmplitude;
curve.figureId = 'figS2a';
curve.featureKind = 'amplitude';
curve.targetValue = targetAmplitude;
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
ampA = pathA.obsAmplitude(idxA);
ampB = pathB.obsAmplitude(idxB);
maxA = pathA.obsMax(idxA);
maxB = pathB.obsMax(idxB);
condA = pathA.directConditionEstimate(idxA);
condB = pathB.directConditionEstimate(idxB);

if ~isempty(matrixA) && ~isempty(matrixB) && norm(matrixA(end, :) - matrixB(1, :), inf) < 1e-10
    matrixB(1, :) = [];
    periodB(1) = [];
    ampB(1) = [];
    maxB(1) = [];
    condB(1) = [];
end

curve = struct();
curve.paramNames = pathA.paramNames;
curve.paramMatrix = [matrixA; matrixB];
curve.Kd = curve.paramMatrix(:, iKd);
curve.AT = curve.paramMatrix(:, iAT);
curve.period = [periodA; periodB];
curve.obsAmplitude = [ampA; ampB];
curve.obsMax = [maxA; maxB];
curve.directConditionEstimate = [condA; condB];
curve.stopReason = {pathA.stopReason, pathB.stopReason};
curve.segments = { ...
    struct('name', 'left', 'Kd', matrixA(:, iKd), 'AT', matrixA(:, iAT), ...
        'period', periodA, 'obsAmplitude', ampA, 'obsMax', maxA, 'directConditionEstimate', condA), ...
    struct('name', 'right', 'Kd', matrixB(:, iKd), 'AT', matrixB(:, iAT), ...
        'period', periodB, 'obsAmplitude', ampB, 'obsMax', maxB, 'directConditionEstimate', condB)};
end
