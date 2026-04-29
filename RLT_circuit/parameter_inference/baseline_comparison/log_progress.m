function log_progress(result, config, evalIdx, elapsedTime)
if ~isfield(config, 'progress') || ~isfield(config.progress, 'enabled') || ...
        ~config.progress.enabled
    return
end

interval = config.progress.interval;
if isempty(interval) || interval < 1
    interval = max(1, floor(config.evalBudget / 20));
end

isFirst = evalIdx == 1;
isInterval = mod(evalIdx, interval) == 0;
isFinal = evalIdx >= config.evalBudget;
if ~(isFirst || isInterval || isFinal)
    return
end

successCount = nnz(result.successFlags(1:evalIdx));
fprintf('[%s] seed=%d eval=%d/%d success=%d bestLoss=%.6e elapsed=%.1fs\n', ...
    result.methodName, result.seed, evalIdx, config.evalBudget, ...
    successCount, result.bestLoss, elapsedTime);
end
