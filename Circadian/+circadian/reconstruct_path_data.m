function path = reconstruct_path_data(model, startParams, seedOrbit, logs, status, ~)
if isempty(logs)
    nPoints = 1;
else
    nPoints = numel(logs) + 1;
end

paramMatrix = NaN(nPoints, numel(startParams));
period = NaN(nPoints, 1);
obsAmplitude = NaN(nPoints, 1);
obsMax = NaN(nPoints, 1);
obsMin = NaN(nPoints, 1);
lambda = NaN(nPoints, 1);
directConditionEstimate = NaN(nPoints, 1);
finalConditionEstimate = NaN(nPoints, 1);
success = true(nPoints, 1);

paramMatrix(1, :) = reshape(startParams, 1, []);
period(1) = seedOrbit.period;
obsAmplitude(1) = seedOrbit.amplitude;
obsMax(1) = seedOrbit.yMax;
obsMin(1) = seedOrbit.yMin;
lambda(1) = 0;

for i = 1:numel(logs)
    paramMatrix(i + 1, :) = reshape(logs(i).params, 1, []);
    period(i + 1) = logs(i).period;
    if isfield(logs, 'obsAmplitude')
        obsAmplitude(i + 1) = logs(i).obsAmplitude;
        obsMax(i + 1) = logs(i).obsMax;
        obsMin(i + 1) = logs(i).obsMin;
    else
        obsAmplitude(i + 1) = logs(i).amplitude;
        obsMax(i + 1) = logs(i).yMax;
        obsMin(i + 1) = logs(i).yMin;
    end
    lambda(i + 1) = logs(i).lambda;
    if isfield(logs, 'directConditionEstimate')
        directConditionEstimate(i + 1) = logs(i).directConditionEstimate;
    end
    if isfield(logs, 'finalConditionEstimate')
        finalConditionEstimate(i + 1) = logs(i).finalConditionEstimate;
    end
    success(i + 1) = all(isfinite([period(i + 1), obsAmplitude(i + 1), obsMax(i + 1), obsMin(i + 1)]));
end

path = struct();
path.paramNames = model.paramNames;
path.paramMatrix = paramMatrix;
path.period = period;
path.obsAmplitude = obsAmplitude;
path.obsMax = obsMax;
path.obsMin = obsMin;
path.lambda = lambda;
path.directConditionEstimate = directConditionEstimate;
path.finalConditionEstimate = finalConditionEstimate;
path.success = success;
path.stopReason = status.reason;
path.stopLambda = status.lambda;
path.stopTriggerValue = status.triggerValue;
end
