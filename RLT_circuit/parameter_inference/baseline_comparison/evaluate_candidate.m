function evalOut = evaluate_candidate(activeParams, config)
activeParams = reshape(activeParams, 1, []);
Parameters = config.baseParameters;
Parameters(config.activeIndex) = activeParams;

evalOut = struct();
evalOut.activeParams = activeParams;
evalOut.Parameters = Parameters;
evalOut.loss = config.penaltyLoss;
evalOut.success = false;
evalOut.msg = '';
evalOut.forwardResult = [];

try
    forwardResult = circuit_forward_orbit(Parameters, config.y0, config.forwardOptions);
    evalOut.forwardResult = forwardResult;

    if ~forwardResult.success
        evalOut.msg = char(forwardResult.msg);
        return
    end

    loss = loss_function(forwardResult.orbit, config.targetOrbit, ...
        config.lossOptions, forwardResult.features, config.targetFeatures);
    if ~isfinite(loss)
        evalOut.msg = 'Non-finite loss.';
        return
    end

    evalOut.loss = loss;
    evalOut.success = true;
    evalOut.msg = '';
catch ME
    evalOut.loss = config.penaltyLoss;
    evalOut.success = false;
    evalOut.msg = ME.message;
end
end
