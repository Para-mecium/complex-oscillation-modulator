function result = baseline_lhs(config, seed)
rng(seed, 'twister');

result = empty_result('lhs', config, seed);
lowerBound = config.activeLowerBound;
upperBound = config.activeUpperBound;
span = upperBound - lowerBound;
numSamples = config.evalBudget;
activeDim = numel(lowerBound);

unitSamples = zeros(numSamples, activeDim);
for dimIdx = 1:activeDim
    order = randperm(numSamples).';
    unitSamples(:, dimIdx) = (order - rand(numSamples, 1)) / numSamples;
end

runTimer = tic;
for evalIdx = 1:numSamples
    activeParams = lowerBound + span .* unitSamples(evalIdx, :);
    evalOut = evaluate_candidate(activeParams, config);
    elapsedTime = toc(runTimer);
    result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
    log_progress(result, config, evalIdx, elapsedTime);
end
result.runtime = toc(runTimer);
result = trim_result(result);
end
