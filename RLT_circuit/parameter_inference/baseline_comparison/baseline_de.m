function result = baseline_de(config, seed)
rng(seed, 'twister');

result = empty_result('de', config, seed);
lowerBound = config.activeLowerBound;
upperBound = config.activeUpperBound;
span = upperBound - lowerBound;
activeDim = numel(lowerBound);

populationSize = config.de.populationSize;
if isempty(populationSize)
    populationSize = min(10 * activeDim, config.evalBudget);
end
populationSize = max(4, min(populationSize, config.evalBudget));

F = config.de.F;
CR = config.de.CR;

population = lowerBound + span .* rand(populationSize, activeDim);
populationLoss = Inf(populationSize, 1);
evalIdx = 0;

runTimer = tic;
for i = 1:populationSize
    evalIdx = evalIdx + 1;
    evalOut = evaluate_candidate(population(i, :), config);
    populationLoss(i) = evalOut.loss;
    elapsedTime = toc(runTimer);
    result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
    log_progress(result, config, evalIdx, elapsedTime);
end

while evalIdx < config.evalBudget
    for targetIdx = 1:populationSize
        if evalIdx >= config.evalBudget
            break
        end

        candidatePool = setdiff(1:populationSize, targetIdx);
        chosen = candidatePool(randperm(numel(candidatePool), 3));
        donor = population(chosen(1), :) + ...
            F * (population(chosen(2), :) - population(chosen(3), :));
        donor = min(max(donor, lowerBound), upperBound);

        trial = population(targetIdx, :);
        crossoverMask = rand(1, activeDim) < CR;
        crossoverMask(randi(activeDim)) = true;
        trial(crossoverMask) = donor(crossoverMask);

        evalIdx = evalIdx + 1;
        evalOut = evaluate_candidate(trial, config);
        elapsedTime = toc(runTimer);
        result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
        log_progress(result, config, evalIdx, elapsedTime);

        if evalOut.loss < populationLoss(targetIdx)
            population(targetIdx, :) = trial;
            populationLoss(targetIdx) = evalOut.loss;
        end
    end
end

result.runtime = toc(runTimer);
result.populationSize = populationSize;
result.F = F;
result.CR = CR;
result.finalPopulation = population;
result.finalPopulationLoss = populationLoss;
result = trim_result(result);
end
