function path = reconstruct_path_data(model, startParams, seedOrbit, logs, status, cfg)
if isempty(logs)
    nPoints = 1;
else
    nPoints = numel(logs) + 1;
end

paramMatrix = NaN(nPoints, numel(startParams));
period = NaN(nPoints, 1);
amplitude = NaN(nPoints, 1);
yMax = NaN(nPoints, 1);
yMin = NaN(nPoints, 1);
lambda = NaN(nPoints, 1);
directConditionEstimate = NaN(nPoints, 1);
finalConditionEstimate = NaN(nPoints, 1);
success = true(nPoints, 1);

paramMatrix(1, :) = reshape(startParams, 1, []);
period(1) = seedOrbit.period;
amplitude(1) = seedOrbit.amplitude;
yMax(1) = seedOrbit.yMax;
yMin(1) = seedOrbit.yMin;
lambda(1) = 0;
for i = 1:numel(logs)
    paramMatrix(i + 1, :) = reshape(logs(i).params, 1, []);
    if isfield(logs, 'period')
        period(i + 1) = logs(i).period;
    end
    if isfield(logs, 'amplitude')
        amplitude(i + 1) = logs(i).amplitude;
    end
    if isfield(logs, 'yMax')
        yMax(i + 1) = logs(i).yMax;
    end
    if isfield(logs, 'yMin')
        yMin(i + 1) = logs(i).yMin;
    end
    lambda(i + 1) = logs(i).lambda;
    if isfield(logs, 'directConditionEstimate')
        directConditionEstimate(i + 1) = logs(i).directConditionEstimate;
    end
    if isfield(logs, 'finalConditionEstimate')
        finalConditionEstimate(i + 1) = logs(i).finalConditionEstimate;
    end
    success(i + 1) = all(isfinite([period(i + 1), amplitude(i + 1), yMax(i + 1), yMin(i + 1)]));
end

path = struct();
path.paramNames = model.paramNames;
path.paramMatrix = paramMatrix;
path.period = period;
path.amplitude = amplitude;
path.yMax = yMax;
path.yMin = yMin;
path.lambda = lambda;
path.directConditionEstimate = directConditionEstimate;
path.finalConditionEstimate = finalConditionEstimate;
path.success = success;
path.stopReason = status.reason;
path.stopLambda = status.lambda;
path.stopTriggerValue = status.triggerValue;
end
