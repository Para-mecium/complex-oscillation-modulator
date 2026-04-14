function cache = build_longevity_modulation_cache(figureId, controlledIdx, controlledNames, targetMinSH, cacheFile)
scriptDir = fileparts(mfilename('fullpath'));
fmamDir = fileparts(scriptDir);

startDir = pwd;
cleanupObj = onCleanup(@() cd(startDir));

obs = [];
M = 75;
PV = struct('name', 'var', 'idx', 1);
params = [30.5, 183, 0.1, 0.1, 3.7, 3.7, 0.3, 0.2, 3.8, 326, 185, 3, 4.8];
y0 = [23.3896986043691; 222.452142059269; 207.783959243676; 216.598138352976];

cd(scriptDir);
System

cd(fmamDir);
addpath(fullfile(fmamDir, 'PO_extract'));
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter); ...
    sys{4}(y(:).', parameter)];
poResult = extract_periodic_orbit(odeFunc, y0, params, struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 400], ...
    'event', 1, ...
    'minCrossings', 3, ...
    'transientFraction', 0, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9)));
if ~poResult.has_orbit
    error('build_longevity_modulation_cache:BaselinePeriodicOrbitNotFound', ...
        'Baseline periodic-orbit extraction failed for %s.', figureId);
end

baselineState = state(obs, params, poResult.orbit_t, poResult.orbit_y, M, PV);
baselineState.updatePeriod();
baselineState.updateVar2();

solverResult = solve_target_modulation( ...
    sys, obs, derivatives, baselineState, controlledIdx, targetMinSH, figureId);
targetView = solverResult.solverView;
targetDerived = solverResult.derivedView;

if isempty(solverResult.logs)
    error('build_longevity_modulation_cache:MissingContinuationLogs', ...
        'No continuation logs were captured for %s.', figureId);
end

cache = struct();
cache.meta = struct();
cache.meta.figureId = char(figureId);
cache.meta.controlledIdx = reshape(double(controlledIdx), 1, []);
cache.meta.controlledNames = reshape(cellstr(controlledNames), 1, []);
cache.meta.targetMinSH = reshape(double(targetMinSH), 1, []);
cache.meta.agingZone = struct('S', 200, 'H', 100);

cache.baseline = export_state_snapshot(baselineState);
cache.target = export_state_snapshot(targetDerived, targetView.params);
cache.path = reconstruct_path_cache(cache.baseline, cache.target, solverResult.logs);

save(cacheFile, '-struct', 'cache');
end

function result = solve_target_modulation(sys, obs, derivatives, baselineState, controlledIdx, targetMinSH, figureId)
targetState = state( ...
    obs, baselineState.params, baselineState.t, baselineState.TS_var, ...
    baselineState.truncationOrder, baselineState.PV);
targetState.updatePeriod();
targetState.updateVar2();
targetView = fmam_state_ops.solverViewFromState(targetState);
targetDerived = fmam_state_ops.derivedViewFromState(targetState);
targetInput = targetView;

lambdaCap = 5e-3;
minLambdaStep = 1e-6;
maxRetries = 16;
acceptedLogs = struct([]);
lambdaOffset = 0;
previousDistance = inf;
finalStatus = struct('completed', false, 'lambda', 0, 'reason', '', 'triggerValue', NaN);
targetReached = false;

for retryIdx = 1:maxRetries
    itemsPer = build_target_items(targetMinSH);

    continuationOptions = struct( ...
        'initialLambdaStep', lambdaCap, ...
        'maxLambdaStep', lambdaCap, ...
        'minLambdaStep', minLambdaStep);
    newtonOptions = struct('maxIterations', 100);

    task = FMAM_ODE(sys, obs, targetInput, itemsPer, controlledIdx, [], 1e-6, ...
        'derivatives', derivatives, ...
        'continuationOptions', continuationOptions, ...
        'newtonOptions', newtonOptions);
    task.isPsiUpdated = true;
    task.needLog = true;
    task.verbose = true;

    task.fit();
    task.step();
    task.errBound = 1e-12;
    task.fit();

    targetView = task.exportSolverView();
    targetDerived = task.exportDerivedView();
    targetInput = targetView;

    acceptedLogs = merge_segment_logs(acceptedLogs, task.logs, lambdaOffset);
    finalStatus = task.continuationStatus;

    currentMinSH = double(targetDerived.varMin(3:4));
    distanceToTarget = norm(currentMinSH - targetMinSH, inf);
    if finalStatus.completed || distanceToTarget <= 1e-6
        targetReached = true;
        break
    end

    segmentLambda = double(finalStatus.lambda);
    if ~(isfinite(segmentLambda) && segmentLambda > 1e-6)
        break
    end

    lambdaOffset = lambdaOffset + (1 - lambdaOffset) * segmentLambda;
    if finalStatus.completed || distanceToTarget >= previousDistance - 1e-6
        lambdaCap = max(lambdaCap * 0.5, 5e-4);
    end
    previousDistance = distanceToTarget;
end

if ~targetReached
    error('build_longevity_modulation_cache:TargetNotReached', ...
        ['Failed to reach deterministic target for %s. Final min(S,H) = ' ...
         '(%.6f, %.6f), status=%s at lambda=%.6f.'], ...
        char(figureId), targetDerived.varMin(3), targetDerived.varMin(4), ...
        char(finalStatus.reason), finalStatus.lambda);
end

odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter); ...
    sys{4}(y(:).', parameter)];
poResult = extract_periodic_orbit(odeFunc, targetDerived.TS_var(1, :).', targetView.params, struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, max([400, 10 * targetDerived.period, targetDerived.t(end)])], ...
    'event', 1, ...
    'minCrossings', 3, ...
    'transientFraction', 0, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9)));
if ~poResult.has_orbit
    error('build_longevity_modulation_cache:TargetPeriodicOrbitNotFound', ...
        'Target periodic-orbit extraction failed for %s.', figureId);
end

targetState = state( ...
    obs, targetView.params, poResult.orbit_t, poResult.orbit_y, ...
    targetInput.truncationOrder, targetView.PV);
targetState.updatePeriod();
targetState.updateVar2();

result = struct( ...
    'state', targetState, ...
    'solverView', fmam_state_ops.solverViewFromState(targetState), ...
    'derivedView', fmam_state_ops.derivedViewFromState(targetState), ...
    'logs', acceptedLogs, ...
    'status', finalStatus);
end

function itemsPer = build_target_items(targetMinSH)
itemsPer = struct([]);
itemsPer(1).prop = 'varMin';
itemsPer(1).idx = 3;
itemsPer(1).target = targetMinSH(1);
itemsPer(2).prop = 'varMin';
itemsPer(2).idx = 4;
itemsPer(2).target = targetMinSH(2);
end

function mergedLogs = merge_segment_logs(baseLogs, segmentLogs, lambdaOffset)
if isempty(segmentLogs)
    mergedLogs = baseLogs;
    return
end

segmentLogs = reshape(segmentLogs, 1, []);
for i = 1:numel(segmentLogs)
    if isfield(segmentLogs, 'lambda') && ~isempty(segmentLogs(i).lambda)
        segmentLogs(i).lambda = lambdaOffset + (1 - lambdaOffset) * double(segmentLogs(i).lambda);
    end
end

if isempty(baseLogs)
    mergedLogs = segmentLogs;
else
    mergedLogs = concat_log_structs(baseLogs, segmentLogs);
end
end

function mergedLogs = concat_log_structs(baseLogs, retryLogs)
fieldNames = union(fieldnames(baseLogs), fieldnames(retryLogs), 'stable');
baseLogs = fill_missing_log_fields(baseLogs, fieldNames);
retryLogs = fill_missing_log_fields(retryLogs, fieldNames);
mergedLogs = [baseLogs, retryLogs];
end

function logs = fill_missing_log_fields(logs, fieldNames)
for i = 1:numel(fieldNames)
    name = fieldNames{i};
    if ~isfield(logs, name)
        [logs.(name)] = deal([]);
    end
end
end

function snapshot = export_state_snapshot(stateLike, params)
snapshot = struct();
if nargin < 2 || isempty(params)
    params = stateLike.params;
end
snapshot.Parameters = reshape(double(params), 1, []);
snapshot.period = double(stateLike.period);
snapshot.varMin = reshape(double(stateLike.varMin), 1, []);
snapshot.varMax = reshape(double(stateLike.varMax), 1, []);
snapshot.t = double(stateLike.t);
snapshot.TS_var = double(stateLike.TS_var);
end

function path = reconstruct_path_cache(baseline, target, logs)
nPoint = numel(logs) + 1;

path = struct();
path.lambda = NaN(nPoint, 1);
path.paramMatrix = NaN(nPoint, numel(baseline.Parameters));
path.minS = NaN(nPoint, 1);
path.minH = NaN(nPoint, 1);

path.lambda(1) = 0;
path.paramMatrix(1, :) = baseline.Parameters;
path.minS(1) = baseline.varMin(3);
path.minH(1) = baseline.varMin(4);

for i = 1:numel(logs)
    path.paramMatrix(i + 1, :) = reshape(double(logs(i).params), 1, []);
    path.minS(i + 1) = extract_logged_property(logs(i), 'varMin_idx_3');
    path.minH(i + 1) = extract_logged_property(logs(i), 'varMin_idx_4');
    if isfield(logs, 'lambda') && ~isempty(logs(i).lambda)
        path.lambda(i + 1) = double(logs(i).lambda);
    end
end

needsTargetAppend = any(abs(path.paramMatrix(end, :) - target.Parameters) > 1e-12) || ...
    abs(path.minS(end) - target.varMin(3)) > 1e-10 || ...
    abs(path.minH(end) - target.varMin(4)) > 1e-10;

if needsTargetAppend
    path.lambda(end + 1, 1) = 1;
    path.paramMatrix(end + 1, :) = target.Parameters;
    path.minS(end + 1, 1) = target.varMin(3);
    path.minH(end + 1, 1) = target.varMin(4);
end

if any(~isfinite(path.lambda))
    path.lambda = estimate_property_progress(path.minS, path.minH, baseline.varMin(3:4), target.varMin(3:4));
end
end

function value = extract_logged_property(logEntry, fieldName)
acceptedField = ['accepted_' fieldName];
if isfield(logEntry, acceptedField) && ~isempty(logEntry.(acceptedField))
    value = double(logEntry.(acceptedField));
else
    value = double(logEntry.(fieldName));
end
end

function lambda = estimate_property_progress(minS, minH, baselineMinSH, targetMinSH)
points = [double(minS(:)), double(minH(:))];
origin = reshape(double(baselineMinSH), 1, []);
direction = reshape(double(targetMinSH), 1, []) - origin;

if norm(direction) <= eps
    lambda = zeros(size(points, 1), 1);
    return
end

rawLambda = ((points - origin) * direction.') / sum(direction.^2);
rawLambda = min(max(rawLambda, 0), 1);
lambda = cummax(rawLambda);
end
