function result = trim_result(result)
n = result.numEvaluations;
result.sampledActiveParams = result.sampledActiveParams(1:n, :);
result.sampledParams = result.sampledParams(1:n, :);
result.losses = result.losses(1:n);
result.bestSoFarLosses = result.bestSoFarLosses(1:n);
result.timeSoFar = result.timeSoFar(1:n);
result.successFlags = result.successFlags(1:n);
result.messages = result.messages(1:n);
end
