function result = baseline_pso(config, seed)
rng(seed, 'twister');

result = empty_result('pso', config, seed);
lowerBound = config.activeLowerBound;
upperBound = config.activeUpperBound;
span = upperBound - lowerBound;
activeDim = numel(lowerBound);

swarmSize = config.pso.swarmSize;
if isempty(swarmSize)
    swarmSize = min(10 * activeDim, config.evalBudget);
end
swarmSize = max(1, min(swarmSize, config.evalBudget));

w = config.pso.inertia;
c1 = config.pso.cognitiveWeight;
c2 = config.pso.socialWeight;
velocityScale = config.pso.initialVelocityScale;

positions = lowerBound + span .* rand(swarmSize, activeDim);
velocities = velocityScale * span .* randn(swarmSize, activeDim);
personalBestPositions = positions;
personalBestLoss = Inf(swarmSize, 1);
globalBestPosition = NaN(1, activeDim);
globalBestLoss = Inf;
evalIdx = 0;

runTimer = tic;
for particleIdx = 1:swarmSize
    evalIdx = evalIdx + 1;
    evalOut = evaluate_candidate(positions(particleIdx, :), config);
    elapsedTime = toc(runTimer);
    result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
    log_progress(result, config, evalIdx, elapsedTime);

    personalBestLoss(particleIdx) = evalOut.loss;
    if evalOut.loss < globalBestLoss
        globalBestLoss = evalOut.loss;
        globalBestPosition = positions(particleIdx, :);
    end
end

while evalIdx < config.evalBudget
    for particleIdx = 1:swarmSize
        if evalIdx >= config.evalBudget
            break
        end

        r1 = rand(1, activeDim);
        r2 = rand(1, activeDim);
        velocities(particleIdx, :) = ...
            w * velocities(particleIdx, :) + ...
            c1 * r1 .* (personalBestPositions(particleIdx, :) - positions(particleIdx, :)) + ...
            c2 * r2 .* (globalBestPosition - positions(particleIdx, :));

        positions(particleIdx, :) = positions(particleIdx, :) + velocities(particleIdx, :);
        below = positions(particleIdx, :) < lowerBound;
        above = positions(particleIdx, :) > upperBound;
        positions(particleIdx, :) = min(max(positions(particleIdx, :), lowerBound), upperBound);
        velocities(particleIdx, below | above) = 0;

        evalIdx = evalIdx + 1;
        evalOut = evaluate_candidate(positions(particleIdx, :), config);
        elapsedTime = toc(runTimer);
        result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
        log_progress(result, config, evalIdx, elapsedTime);

        if evalOut.loss < personalBestLoss(particleIdx)
            personalBestLoss(particleIdx) = evalOut.loss;
            personalBestPositions(particleIdx, :) = positions(particleIdx, :);
        end
        if evalOut.loss < globalBestLoss
            globalBestLoss = evalOut.loss;
            globalBestPosition = positions(particleIdx, :);
        end
    end
end

result.runtime = toc(runTimer);
result.swarmSize = swarmSize;
result.inertia = w;
result.cognitiveWeight = c1;
result.socialWeight = c2;
result.finalPositions = positions;
result.finalVelocities = velocities;
result.personalBestPositions = personalBestPositions;
result.personalBestLoss = personalBestLoss;
result = trim_result(result);
end
