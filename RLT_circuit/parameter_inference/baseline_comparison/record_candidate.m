function result = record_candidate(result, evalIdx, evalOut, config, elapsedTime)
if nargin < 5
    elapsedTime = NaN;
end

result.sampledActiveParams(evalIdx, :) = evalOut.activeParams;
result.sampledParams(evalIdx, :) = evalOut.Parameters;
result.losses(evalIdx) = evalOut.loss;
result.successFlags(evalIdx) = evalOut.success;
result.messages{evalIdx} = evalOut.msg;
result.timeSoFar(evalIdx) = elapsedTime;
result.numEvaluations = evalIdx;

if evalOut.loss < result.bestLoss
    result.bestLoss = evalOut.loss;
    result.bestEvalIndex = evalIdx;
    result.bestActiveParams = evalOut.activeParams;
    result.bestParameters = evalOut.Parameters;

    if config.saveBestOrbit
        result.bestForwardResult = evalOut.forwardResult;
    end
end

result.bestSoFarLosses(evalIdx) = result.bestLoss;
end
