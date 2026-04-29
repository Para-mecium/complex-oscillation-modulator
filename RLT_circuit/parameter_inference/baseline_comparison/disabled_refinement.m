function refinement = disabled_refinement(config)
activeDim = numel(config.activeIndex);
paramDim = numel(config.baseParameters);

refinement = struct();
refinement.enabled = false;
refinement.method = 'none';
refinement.budget = 0;
refinement.topK = 0;
refinement.numEvaluations = 0;
refinement.runtime = 0;
refinement.startEvalIndices = [];
refinement.startActiveParams = NaN(0, activeDim);
refinement.startLosses = NaN(0, 1);
refinement.bestLoss = NaN;
refinement.bestActiveParams = NaN(1, activeDim);
refinement.bestParameters = NaN(1, paramDim);
refinement.bestForwardResult = [];
refinement.traceActiveParams = NaN(0, activeDim);
refinement.traceParams = NaN(0, paramDim);
refinement.traceLosses = NaN(0, 1);
refinement.bestSoFarLosses = NaN(0, 1);
refinement.globalBestSoFarLosses = NaN(0, 1);
refinement.timeSoFar = NaN(0, 1);
refinement.successFlags = false(0, 1);
refinement.messages = cell(0, 1);
refinement.startSummaries = struct([]);
end
