function result = empty_result(methodName, config, seed)
activeDim = numel(config.activeIndex);
paramDim = numel(config.baseParameters);
budget = config.evalBudget;

result = struct();
result.methodName = methodName;
result.seed = seed;
result.evalBudget = budget;
result.numEvaluations = 0;
result.runtime = NaN;

result.activeIndex = config.activeIndex;
result.activeLowerBound = config.activeLowerBound;
result.activeUpperBound = config.activeUpperBound;
result.activeNames = config.activeNames;

result.sampledActiveParams = NaN(budget, activeDim);
result.sampledParams = NaN(budget, paramDim);
result.losses = NaN(budget, 1);
result.bestSoFarLosses = NaN(budget, 1);
result.timeSoFar = NaN(budget, 1);
result.successFlags = false(budget, 1);
result.messages = cell(budget, 1);

result.bestLoss = Inf;
result.bestEvalIndex = NaN;
result.bestActiveParams = NaN(1, activeDim);
result.bestParameters = NaN(1, paramDim);
result.bestForwardResult = [];
end
