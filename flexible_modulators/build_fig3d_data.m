function data = build_fig3d_data(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

amplitudes = cfg.fig3d.amplitudes;
periods = cfg.fig3d.periods;
assert(numel(amplitudes) == numel(periods), ...
    'cfg.fig3d.amplitudes and cfg.fig3d.periods must align.');

curves = cell(1, numel(amplitudes));
curveFiles = cell(1, numel(amplitudes));
for i = 1:numel(amplitudes)
    curveFiles{i} = flexmod.curve_cache_file(cfg, 'fig3d', 'amplitude', amplitudes(i));
    curves{i} = flexmod.cache_get_or_create(curveFiles{i}, ...
        @() build_iso_amplitude_curve(cfg, amplitudes(i), periods(i)));
end

data = struct();
data.curves = curves;
data.curveFiles = curveFiles;
end

function curve = build_iso_amplitude_curve(cfg, targetAmplitude, targetPeriod)
seedCfg = struct( ...
    'modelType', 'base', ...
    'startParams', cfg.base.initialParams, ...
    'goalOrder', {{'amplitude', 'period'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'period', targetPeriod), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', 1));
seed = run_modulation_task(seedCfg);

left = run_modulation_task(struct( ...
    'modelType', 'base', ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'amplitude', 'I1'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'I1', cfg.fig3d.I1Range(1)), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));
left = retry_with_et_continuation(cfg, left, targetAmplitude);

right = run_modulation_task(struct( ...
    'modelType', 'base', ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'amplitude', 'I1'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'I1', cfg.fig3d.I1Range(2)), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));
right = retry_with_et_continuation(cfg, right, targetAmplitude);

leftPath = flexmod.trim_path(left.path, 'iso_amplitude', cfg);
rightPath = flexmod.trim_path(right.path, 'iso_amplitude', cfg);
curve = combine_branches(leftPath, rightPath, targetAmplitude);
curve.figureId = 'fig3d';
curve.featureKind = 'amplitude';
curve.targetValue = targetAmplitude;
curve.targetPeriodSeed = targetPeriod;
end

function result = retry_with_et_continuation(cfg, result, targetAmplitude)
if strcmp(result.path.stopReason, 'target_reached')
    return
end

[~, ok] = estimate_retry_et_target(result.path);
if ~ok
    return
end

retryStepCap = max(cfg.fmam.dlambdaCap * 0.25, 1e-3);
retry = run_modulation_task(struct( ...
    'modelType', 'base', ...
    'initialSolverView', result.solverView, ...
    'initialDerivedView', result.derivedView, ...
    'startParams', result.params, ...
    'goalOrder', {{'amplitude', 'I1'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'I1', cfg.fig3d.I1Range(1)), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', retryStepCap)));

result = merge_retry_result(result, retry);
end

function [targetET, ok] = estimate_retry_et_target(path)
targetET = NaN;
ok = false;

iET = flexmod.parameter_index(path, 'ET');
valid = find(all(isfinite(path.paramMatrix), 2) & isfinite(path.period) & isfinite(path.amplitude));
if numel(valid) < 2
    return
end

tailCount = min(6, numel(valid));
tailIdx = valid(end - tailCount + 1:end);
tailET = path.paramMatrix(tailIdx, iET);
dET = median(diff(tailET));
if ~(isfinite(dET) && abs(dET) > 0)
    dET = (tailET(end) - tailET(1)) / max(tailCount - 1, 1);
end
if ~(isfinite(dET) && abs(dET) > 0)
    return
end

etIncrement = min(max(4 * abs(dET), 0.05), 0.25) * sign(dET);
targetET = tailET(end) + etIncrement;
ok = isfinite(targetET);
end

function merged = merge_retry_result(base, retry)
merged = retry;
merged.path = merge_retry_path(base.path, retry.path);
if isfield(base, 'logs') && isfield(retry, 'logs')
    if isempty(base.logs)
        merged.logs = retry.logs;
    elseif isempty(retry.logs)
        merged.logs = base.logs;
    else
        merged.logs = concat_log_structs(base.logs, retry.logs);
    end
end
end

function mergedLogs = concat_log_structs(baseLogs, retryLogs)
fieldNames = union(fieldnames(baseLogs), fieldnames(retryLogs), 'stable');
baseLogs = fill_missing_log_fields(baseLogs, fieldNames);
retryLogs = fill_missing_log_fields(retryLogs, fieldNames);
mergedLogs = [baseLogs retryLogs];
end

function logs = fill_missing_log_fields(logs, fieldNames)
for i = 1:numel(fieldNames)
    name = fieldNames{i};
    if ~isfield(logs, name)
        [logs.(name)] = deal([]);
    end
end
end

function mergedPath = merge_retry_path(basePath, retryPath)
mergedPath = basePath;
fields = {'paramMatrix', 'period', 'amplitude', 'yMax', 'yMin', 'lambda', ...
    'directConditionEstimate', 'finalConditionEstimate', 'success'};

for i = 1:numel(fields)
    name = fields{i};
    retryValues = retryPath.(name);
    if size(retryValues, 1) <= 1
        continue
    end

    if strcmp(name, 'lambda')
        lambda0 = basePath.lambda(end);
        retryValues = lambda0 + (1 - lambda0) * retryValues;
    end

    mergedPath.(name) = [basePath.(name); retryValues(2:end, :)];
end

mergedPath.stopReason = retryPath.stopReason;
mergedPath.stopLambda = retryPath.stopLambda;
mergedPath.stopTriggerValue = retryPath.stopTriggerValue;
end

function curve = combine_branches(leftPath, rightPath, targetAmplitude)
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
curve.targetAmplitude = targetAmplitude;
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
