function data = build_fig3c_data(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

periods = cfg.fig3c.periods;
amps = cfg.fig3c.amplitudes;
assert(numel(periods) == numel(amps), ...
    'cfg.fig3c.periods and cfg.fig3c.amplitudes must align.');

curves = cell(1, numel(periods));
curveFiles = cell(1, numel(periods));
for i = 1:numel(periods)
    curveFiles{i} = flexmod.curve_cache_file(cfg, 'fig3c', 'period', periods(i));
    curves{i} = flexmod.cache_get_or_create(curveFiles{i}, ...
        @() build_iso_period_curve(cfg, periods(i), amps(i)));
end

data = struct();
data.curves = curves;
data.curveFiles = curveFiles;
end

function curve = build_iso_period_curve(cfg, targetPeriod, targetAmp)
seedCfg = struct( ...
    'modelType', 'base', ...
    'startParams', cfg.base.initialParams, ...
    'goalOrder', {{'period', 'amplitude'}}, ...
    'goals', struct('period', targetPeriod, 'amplitude', targetAmp), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', 1));
seed = run_modulation_task(seedCfg);

left = run_modulation_task(struct( ...
    'modelType', 'base', ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'period', 'I1'}}, ...
    'goals', struct('period', targetPeriod, 'I1', cfg.fig3c.I1Range(1)), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

right = run_modulation_task(struct( ...
    'modelType', 'base', ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'period', 'I1'}}, ...
    'goals', struct('period', targetPeriod, 'I1', cfg.fig3c.I1Range(2)), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

leftPath = flexmod.trim_path(left.path, 'iso_period', cfg);
rightPath = flexmod.trim_path(right.path, 'iso_period', cfg);
curve = combine_branches(leftPath, rightPath, targetPeriod);
curve.figureId = 'fig3c';
curve.featureKind = 'period';
curve.targetValue = targetPeriod;
curve.targetAmplitudeSeed = targetAmp;
end

function curve = combine_branches(leftPath, rightPath, targetPeriod)
iI1 = flexmod.parameter_index(leftPath, 'I1');
iET = flexmod.parameter_index(leftPath, 'ET');

leftOrder = sortrows([(1:size(leftPath.paramMatrix, 1)).', leftPath.paramMatrix(:, iI1)], 2);
rightOrder = sortrows([(1:size(rightPath.paramMatrix, 1)).', rightPath.paramMatrix(:, iI1)], 2);
leftIdx = leftOrder(:, 1);
rightIdx = rightOrder(:, 1);

leftMatrix = leftPath.paramMatrix(leftIdx, :);
leftPeriod = leftPath.period(leftIdx);
leftAmplitude = leftPath.amplitude(leftIdx);
leftCond = leftPath.directConditionEstimate(leftIdx);

rightMatrix = rightPath.paramMatrix(rightIdx, :);
rightPeriod = rightPath.period(rightIdx);
rightAmplitude = rightPath.amplitude(rightIdx);
rightCond = rightPath.directConditionEstimate(rightIdx);

if ~isempty(leftMatrix) && ~isempty(rightMatrix) && norm(leftMatrix(end, :) - rightMatrix(1, :), inf) < 1e-10
    rightMatrix(1, :) = [];
    rightPeriod(1) = [];
    rightAmplitude(1) = [];
    rightCond(1) = [];
end

curve = struct();
curve.targetPeriod = targetPeriod;
curve.paramNames = leftPath.paramNames;
curve.paramMatrix = [leftMatrix; rightMatrix];
curve.I1 = curve.paramMatrix(:, iI1);
curve.ET = curve.paramMatrix(:, iET);
curve.period = [leftPeriod; rightPeriod];
curve.amplitude = [leftAmplitude; rightAmplitude];
curve.directConditionEstimate = [leftCond; rightCond];
curve.stopReason = {leftPath.stopReason, rightPath.stopReason};
curve.segments = { ...
    struct('name', 'left', 'I1', leftMatrix(:, iI1), 'ET', leftMatrix(:, iET), ...
        'period', leftPeriod, 'amplitude', leftAmplitude, 'directConditionEstimate', leftCond), ...
    struct('name', 'right', 'I1', rightMatrix(:, iI1), 'ET', rightMatrix(:, iET), ...
        'period', rightPeriod, 'amplitude', rightAmplitude, 'directConditionEstimate', rightCond)};
end
