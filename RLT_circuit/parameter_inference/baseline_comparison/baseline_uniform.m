function result = baseline_uniform(config, seed)
rng(seed, 'twister');

result = empty_result('uniform', config, seed);
lowerBound = config.activeLowerBound;
upperBound = config.activeUpperBound;
span = upperBound - lowerBound;

runTimer = tic;
for evalIdx = 1:config.evalBudget
    activeParams = lowerBound + span .* rand(1, numel(lowerBound));
    evalOut = evaluate_candidate(activeParams, config);
    elapsedTime = toc(runTimer);
    result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
    log_progress(result, config, evalIdx, elapsedTime);
end
result.runtime = toc(runTimer);
result = trim_result(result);
end
