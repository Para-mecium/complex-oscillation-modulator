function refinement = run_fmincon(config, seed, baselineResult)
rng(seed, 'twister');

if exist('fmincon', 'file') ~= 2
    error('run_fmincon:MissingOptimizer', ...
        'fmincon is not available. Install or enable Optimization Toolbox.');
end

activeDim = numel(config.activeIndex);
paramDim = numel(config.baseParameters);
budget = config.refinement.budget;
topK = config.refinement.topK;

refinement = struct();
refinement.enabled = true;
refinement.method = 'fmincon';
refinement.budget = budget;
refinement.topK = topK;
refinement.numEvaluations = 0;
refinement.runtime = NaN;
refinement.startEvalIndices = [];
refinement.startActiveParams = NaN(0, activeDim);
refinement.startLosses = NaN(0, 1);
refinement.bestLoss = baselineResult.bestLoss;
refinement.bestActiveParams = baselineResult.bestActiveParams;
refinement.bestParameters = baselineResult.bestParameters;
refinement.bestForwardResult = [];
refinement.traceActiveParams = NaN(budget, activeDim);
refinement.traceParams = NaN(budget, paramDim);
refinement.traceLosses = NaN(budget, 1);
refinement.bestSoFarLosses = NaN(budget, 1);
refinement.globalBestSoFarLosses = NaN(budget, 1);
refinement.timeSoFar = NaN(budget, 1);
refinement.successFlags = false(budget, 1);
refinement.messages = cell(budget, 1);
refinement.startSummaries = repmat(empty_start_summary(activeDim), max(topK, 0), 1);

if budget <= 0 || topK <= 0
    refinement.runtime = 0;
    refinement.startSummaries = struct([]);
    refinement = trim_refinement(refinement);
    return
end

startIdx = select_start_indices(baselineResult, topK);
numStarts = numel(startIdx);
refinement.startEvalIndices = startIdx(:).';
refinement.startActiveParams = baselineResult.sampledActiveParams(startIdx, :);
refinement.startLosses = baselineResult.losses(startIdx);
refinement.startSummaries = repmat(empty_start_summary(activeDim), numStarts, 1);

fprintf('[%s/refine:fmincon] seed=%d starts: topK=%d budget=%d\n', ...
    baselineResult.methodName, seed, numStarts, budget);

runTimer = tic;
for startOrder = 1:numStarts
    if refinement.numEvaluations >= budget
        break
    end

    startX = baselineResult.sampledActiveParams(startIdx(startOrder), :);
    startLoss = baselineResult.losses(startIdx(startOrder));
    [refinement, startSummary] = refine_one_start( ...
        refinement, config, baselineResult.methodName, startOrder, ...
        startIdx(startOrder), startX, startLoss, runTimer);
    refinement.startSummaries(startOrder) = startSummary;
end

refinement.runtime = toc(runTimer);
refinement = trim_refinement(refinement);
end

function [refinement, startSummary] = refine_one_start( ...
    refinement, config, methodName, startOrder, startEvalIndex, startX, ...
    startLoss, runTimer)

activeDim = numel(startX);
startEvalCount = refinement.numEvaluations;
cacheX = startX;
cacheLoss = startLoss;

startSummary = empty_start_summary(activeDim);
startSummary.startOrder = startOrder;
startSummary.startEvalIndex = startEvalIndex;
startSummary.startActiveParams = startX;
startSummary.startLoss = startLoss;
startSummary.bestActiveParams = startX;
startSummary.bestLoss = startLoss;

remainingBudget = refinement.budget - refinement.numEvaluations;
options = optimoptions('fmincon', ...
    'Algorithm', config.refinement.fminconAlgorithm, ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', max(1, remainingBudget + 1), ...
    'FiniteDifferenceType', 'forward', ...
    'StepTolerance', config.refinement.fminconStepTolerance, ...
    'OptimalityTolerance', config.refinement.fminconOptimalityTolerance, ...
    'OutputFcn', @stop_when_budget_used);

try
    [xOpt, fOpt, exitflag, output] = fmincon( ...
        @objective, startX, [], [], [], [], ...
        config.activeLowerBound, config.activeUpperBound, [], options);
    startSummary.exitflag = exitflag;
    startSummary.output = output;
    if fOpt < startSummary.bestLoss
        startSummary.bestActiveParams = xOpt;
        startSummary.bestLoss = fOpt;
    end
catch ME
    startSummary.exitflag = NaN;
    startSummary.output = struct('message', ME.message);
    startSummary.message = ME.message;
end

startSummary.numEvaluations = refinement.numEvaluations - startEvalCount;

    function loss = objective(x)
        x = reshape(x, 1, []);
        x = min(max(x, config.activeLowerBound), config.activeUpperBound);

        cachedIdx = find_cached_point(cacheX, x);
        if ~isempty(cachedIdx)
            loss = cacheLoss(cachedIdx);
            return
        end

        if refinement.numEvaluations >= refinement.budget
            loss = refinement.bestLoss;
            return
        end

        evalOut = evaluate_candidate(x, config);
        refinement = record_refinement_candidate(refinement, evalOut, config, runTimer);
        loss = evalOut.loss;

        cacheX(end + 1, :) = x;
        cacheLoss(end + 1, 1) = loss;

        if loss < startSummary.bestLoss
            startSummary.bestLoss = loss;
            startSummary.bestActiveParams = x;
        end

        log_refinement_progress(refinement, config, methodName, runTimer);
    end

    function stop = stop_when_budget_used(~, ~, ~)
        stop = refinement.numEvaluations >= refinement.budget;
    end
end

function cachedIdx = find_cached_point(cacheX, x)
if isempty(cacheX)
    cachedIdx = [];
    return
end

tolerance = 1e-12;
distance = max(abs(cacheX - x), [], 2);
cachedIdx = find(distance <= tolerance, 1, 'first');
end

function refinement = record_refinement_candidate(refinement, evalOut, config, runTimer)
evalIdx = refinement.numEvaluations + 1;
refinement.traceActiveParams(evalIdx, :) = evalOut.activeParams;
refinement.traceParams(evalIdx, :) = evalOut.Parameters;
refinement.traceLosses(evalIdx) = evalOut.loss;
refinement.successFlags(evalIdx) = evalOut.success;
refinement.messages{evalIdx} = evalOut.msg;
refinement.timeSoFar(evalIdx) = toc(runTimer);
refinement.numEvaluations = evalIdx;

if evalOut.loss < refinement.bestLoss
    refinement.bestLoss = evalOut.loss;
    refinement.bestActiveParams = evalOut.activeParams;
    refinement.bestParameters = evalOut.Parameters;
    if config.saveBestOrbit
        refinement.bestForwardResult = evalOut.forwardResult;
    end
end

refinement.globalBestSoFarLosses(evalIdx) = refinement.bestLoss;
if evalIdx == 1
    refinement.bestSoFarLosses(evalIdx) = evalOut.loss;
else
    refinement.bestSoFarLosses(evalIdx) = min( ...
        refinement.bestSoFarLosses(evalIdx - 1), evalOut.loss);
end
end

function log_refinement_progress(refinement, config, methodName, runTimer)
if ~isfield(config, 'progress') || ~isfield(config.progress, 'enabled') || ...
        ~config.progress.enabled
    return
end

interval = config.progress.interval;
if isempty(interval) || interval < 1
    interval = max(1, floor(refinement.budget / 20));
end

evalIdx = refinement.numEvaluations;
if ~(evalIdx == 1 || mod(evalIdx, interval) == 0 || evalIdx >= refinement.budget)
    return
end

successCount = nnz(refinement.successFlags(1:evalIdx));
fprintf(['[%s/refine:fmincon] eval=%d/%d success=%d bestLoss=%.6e ' ...
    'elapsed=%.1fs\n'], methodName, evalIdx, refinement.budget, ...
    successCount, refinement.bestLoss, toc(runTimer));
end

function startIdx = select_start_indices(baselineResult, topK)
validMask = baselineResult.successFlags & isfinite(baselineResult.losses);
validIdx = find(validMask);
if isempty(validIdx)
    validIdx = find(isfinite(baselineResult.losses));
end
if isempty(validIdx)
    startIdx = [];
    return
end

[~, order] = sort(baselineResult.losses(validIdx), 'ascend');
sortedIdx = validIdx(order);
keep = true(size(sortedIdx));
for i = 1:numel(sortedIdx)
    if ~keep(i)
        continue
    end
    for j = (i + 1):numel(sortedIdx)
        if isequal(baselineResult.sampledActiveParams(sortedIdx(i), :), ...
                baselineResult.sampledActiveParams(sortedIdx(j), :))
            keep(j) = false;
        end
    end
end
sortedIdx = sortedIdx(keep);
startIdx = sortedIdx(1:min(topK, numel(sortedIdx)));
end

function summary = empty_start_summary(activeDim)
summary = struct();
summary.startOrder = NaN;
summary.startEvalIndex = NaN;
summary.startActiveParams = NaN(1, activeDim);
summary.startLoss = NaN;
summary.numEvaluations = 0;
summary.bestActiveParams = NaN(1, activeDim);
summary.bestLoss = NaN;
summary.exitflag = NaN;
summary.output = struct();
summary.message = '';
end

function refinement = trim_refinement(refinement)
n = refinement.numEvaluations;
refinement.traceActiveParams = refinement.traceActiveParams(1:n, :);
refinement.traceParams = refinement.traceParams(1:n, :);
refinement.traceLosses = refinement.traceLosses(1:n);
refinement.bestSoFarLosses = refinement.bestSoFarLosses(1:n);
refinement.globalBestSoFarLosses = refinement.globalBestSoFarLosses(1:n);
refinement.timeSoFar = refinement.timeSoFar(1:n);
refinement.successFlags = refinement.successFlags(1:n);
refinement.messages = refinement.messages(1:n);
end
