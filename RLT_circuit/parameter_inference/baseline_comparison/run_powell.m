function refinement = run_powell(config, seed, baselineResult)
rng(seed, 'twister');

activeDim = numel(config.activeIndex);
paramDim = numel(config.baseParameters);
budget = config.refinement.budget;
topK = config.refinement.topK;

refinement = struct();
refinement.enabled = true;
refinement.method = 'powell';
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

fprintf('[%s/refine:powell] seed=%d starts: topK=%d budget=%d\n', ...
    baselineResult.methodName, seed, numStarts, budget);

runTimer = tic;
lowerBound = config.activeLowerBound;
upperBound = config.activeUpperBound;
span = upperBound - lowerBound;

for startOrder = 1:numStarts
    if refinement.numEvaluations >= budget
        break
    end

    startX = baselineResult.sampledActiveParams(startIdx(startOrder), :);
    startLoss = baselineResult.losses(startIdx(startOrder));
    [refinement, startSummary] = refine_one_start( ...
        refinement, config, baselineResult.methodName, startOrder, ...
        startIdx(startOrder), startX, startLoss, lowerBound, upperBound, ...
        span, runTimer);
    refinement.startSummaries(startOrder) = startSummary;
end

refinement.runtime = toc(runTimer);
refinement = trim_refinement(refinement);
end

function [refinement, startSummary] = refine_one_start( ...
    refinement, config, methodName, startOrder, startEvalIndex, startX, ...
    startLoss, lowerBound, upperBound, span, runTimer)

activeDim = numel(startX);
directions = diag(span);
x = min(max(startX, lowerBound), upperBound);
xLoss = startLoss;
startEvalCount = refinement.numEvaluations;
maxSweeps = config.refinement.maxSweeps;

startSummary = empty_start_summary(activeDim);
startSummary.startOrder = startOrder;
startSummary.startEvalIndex = startEvalIndex;
startSummary.startActiveParams = startX;
startSummary.startLoss = startLoss;
startSummary.bestActiveParams = x;
startSummary.bestLoss = xLoss;

for sweepIdx = 1:maxSweeps
    if refinement.numEvaluations >= refinement.budget
        break
    end

    sweepStartX = x;
    sweepStartLoss = xLoss;
    largestDecrease = -Inf;
    largestDecreaseIdx = 1;

    for directionIdx = 1:activeDim
        if refinement.numEvaluations >= refinement.budget
            break
        end

        oldLoss = xLoss;
        [x, xLoss, refinement] = line_search( ...
            x, directions(directionIdx, :), xLoss, refinement, config, ...
            methodName, runTimer);
        decrease = oldLoss - xLoss;
        if decrease > largestDecrease
            largestDecrease = decrease;
            largestDecreaseIdx = directionIdx;
        end
    end

    newDirection = x - sweepStartX;
    scaledStep = norm(newDirection ./ span);
    lossDecrease = sweepStartLoss - xLoss;

    if any(newDirection ~= 0) && refinement.numEvaluations < refinement.budget
        [x, xLoss, refinement] = line_search( ...
            x, newDirection, xLoss, refinement, config, methodName, runTimer);
        directions(largestDecreaseIdx, :) = newDirection;
    end

    if scaledStep < config.refinement.tolerance || ...
            abs(lossDecrease) < config.refinement.lossTolerance
        break
    end
end

startSummary.numEvaluations = refinement.numEvaluations - startEvalCount;
startSummary.bestActiveParams = x;
startSummary.bestLoss = xLoss;
end

function [xBest, fBest, refinement] = line_search( ...
    x0, direction, f0, refinement, config, methodName, runTimer)

if norm(direction) == 0
    xBest = x0;
    fBest = f0;
    return
end

[alphaLower, alphaUpper] = bounded_alpha_interval( ...
    x0, direction, config.activeLowerBound, config.activeUpperBound);
if alphaUpper <= alphaLower
    xBest = x0;
    fBest = f0;
    return
end

lineBudget = min(config.refinement.lineSearchMaxEval, ...
    refinement.budget - refinement.numEvaluations);
if lineBudget <= 0
    xBest = x0;
    fBest = f0;
    return
end

golden = (sqrt(5) - 1) / 2;
a = alphaLower;
b = alphaUpper;
c = b - golden * (b - a);
d = a + golden * (b - a);

[fc, refinement] = evaluate_alpha(c, x0, direction, refinement, config, ...
    methodName, runTimer);
if refinement.numEvaluations >= refinement.budget || lineBudget <= 1
    [xBest, fBest] = best_line_point(x0, direction, f0, c, fc);
    return
end

[fd, refinement] = evaluate_alpha(d, x0, direction, refinement, config, ...
    methodName, runTimer);
used = 2;

while used < lineBudget && refinement.numEvaluations < refinement.budget && ...
        abs(b - a) > config.refinement.minLineSearchWidth
    if fc < fd
        b = d;
        d = c;
        fd = fc;
        c = b - golden * (b - a);
        [fc, refinement] = evaluate_alpha(c, x0, direction, refinement, ...
            config, methodName, runTimer);
    else
        a = c;
        c = d;
        fc = fd;
        d = a + golden * (b - a);
        [fd, refinement] = evaluate_alpha(d, x0, direction, refinement, ...
            config, methodName, runTimer);
    end
    used = used + 1;
end

[xBest, fBest] = best_line_point(x0, direction, f0, c, fc, d, fd);
end

function [loss, refinement] = evaluate_alpha( ...
    alpha, x0, direction, refinement, config, methodName, runTimer)

activeParams = x0 + alpha * direction;
activeParams = min(max(activeParams, config.activeLowerBound), ...
    config.activeUpperBound);
evalOut = evaluate_candidate(activeParams, config);
refinement = record_refinement_candidate(refinement, evalOut, config, runTimer);
loss = evalOut.loss;
log_refinement_progress(refinement, config, methodName, runTimer);
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
fprintf(['[%s/refine:powell] eval=%d/%d success=%d bestLoss=%.6e ' ...
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

function [alphaLower, alphaUpper] = bounded_alpha_interval(x, direction, lowerBound, upperBound)
alphaLower = -Inf;
alphaUpper = Inf;
for i = 1:numel(x)
    if direction(i) > 0
        alphaLower = max(alphaLower, (lowerBound(i) - x(i)) / direction(i));
        alphaUpper = min(alphaUpper, (upperBound(i) - x(i)) / direction(i));
    elseif direction(i) < 0
        alphaLower = max(alphaLower, (upperBound(i) - x(i)) / direction(i));
        alphaUpper = min(alphaUpper, (lowerBound(i) - x(i)) / direction(i));
    end
end
end

function [xBest, fBest] = best_line_point(x0, direction, f0, varargin)
fBest = f0;
xBest = x0;
for i = 1:2:numel(varargin)
    alpha = varargin{i};
    loss = varargin{i + 1};
    if loss < fBest
        fBest = loss;
        xBest = x0 + alpha * direction;
    end
end
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
