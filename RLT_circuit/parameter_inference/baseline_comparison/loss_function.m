function loss = loss_function(candidateOrbit, targetOrbit, options, candidateFeatures, targetFeatures)
if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'name') || isempty(options.name)
    options.name = 'relative_l2_orbit';
end
if ~isfield(options, 'compareNumPoints') || isempty(options.compareNumPoints)
    options.compareNumPoints = size(targetOrbit.y, 1);
end
if ~isfield(options, 'phaseAlignment') || isempty(options.phaseAlignment)
    options.phaseAlignment = true;
end
switch lower(options.name)
    case {'relative_l2_orbit', 'relative_l2'}
        loss = loss_relative_l2_orbit(candidateOrbit, targetOrbit, options);
    case 'mse'
        loss = loss_mse(candidateOrbit, targetOrbit, options);
    case 'rmse'
        loss = loss_rmse(candidateOrbit, targetOrbit, options);
    case {'property_difference', 'property_diff'}
        loss = loss_property_difference( ...
            candidateOrbit, targetOrbit, options, candidateFeatures, targetFeatures);
    otherwise
        error('loss_function:UnknownLoss', 'Unknown loss name: %s.', options.name);
end

if ~isfinite(loss)
    loss = Inf;
end
end

function loss = loss_relative_l2_orbit(candidateOrbit, targetOrbit, options)
[candidateY, targetY] = resample_orbits(candidateOrbit, targetOrbit, options);
diffY = candidateY - targetY;

denom = norm(targetY(:));
if denom == 0
    denom = 1;
end
loss = norm(diffY(:)) / denom;
loss = add_period_loss(loss, candidateOrbit, targetOrbit, options);
end

function loss = loss_mse(candidateOrbit, targetOrbit, options)
[candidateY, targetY] = resample_orbits(candidateOrbit, targetOrbit, options);
diffY = candidateY - targetY;
loss = mean(diffY(:).^2);
loss = add_period_loss(loss, candidateOrbit, targetOrbit, options);
end

function loss = loss_rmse(candidateOrbit, targetOrbit, options)
[candidateY, targetY] = resample_orbits(candidateOrbit, targetOrbit, options);
diffY = candidateY - targetY;
loss = sqrt(mean(diffY(:).^2));
loss = add_period_loss(loss, candidateOrbit, targetOrbit, options);
end

function loss = loss_property_difference(candidateOrbit, targetOrbit, options, ...
    candidateFeatures, targetFeatures)
stateIndex = options.propertyStateIndex;
candidateProps = [ ...
    candidateFeatures.period, ...
    reshape(candidateFeatures.state.amplitude(stateIndex), 1, [])];
targetProps = [ ...
    targetFeatures.period, ...
    reshape(targetFeatures.state.amplitude(stateIndex), 1, [])];

scaleFloor = options.propertyScaleFloor;
scale = max(abs(targetProps), scaleFloor);
relativeDiff = (candidateProps - targetProps) ./ scale;

weights = options.propertyWeights;
loss = sqrt(mean((weights .* relativeDiff).^2));
end


%%
function loss = add_period_loss(loss, candidateOrbit, targetOrbit, options)
if ~isfield(options, 'periodWeight') || isempty(options.periodWeight)
    options.periodWeight = 0;
end

if options.periodWeight > 0 && isfield(candidateOrbit, 'period') && isfield(targetOrbit, 'period')
    periodScale = abs(targetOrbit.period);
    if periodScale == 0
        periodScale = 1;
    end
    periodLoss = abs(candidateOrbit.period - targetOrbit.period) / periodScale;
    loss = loss + options.periodWeight * periodLoss;
end
end

function [candidateY, targetY] = resample_orbits(candidateOrbit, targetOrbit, options)
nPoints = options.compareNumPoints;
targetPhase = normalized_phase(targetOrbit.t);
candidatePhase = normalized_phase(candidateOrbit.t);
commonPhase = linspace(0, 1, nPoints + 1).';
commonPhase(end) = [];

targetY = interp1(targetPhase, targetOrbit.y, commonPhase, 'linear');
candidateY = interp1(candidatePhase, candidateOrbit.y, commonPhase, 'linear');

if any(~isfinite(targetY), 'all') || any(~isfinite(candidateY), 'all')
    error('loss_function:NonFiniteOrbit', ...
        'Periodic-orbit samples contain non-finite values.');
end

if options.phaseAlignment
    candidateY = align_candidate_phase(candidateY, targetY);
end
end

function phase = normalized_phase(t)
t = t(:);
if numel(t) < 2 || t(end) == t(1)
    error('loss_function:InvalidOrbitTime', ...
        'Orbit time vector must contain at least two distinct values.');
end
phase = (t - t(1)) / (t(end) - t(1));
end

function alignedCandidateY = align_candidate_phase(candidateY, targetY)
nPoints = size(candidateY, 1);
targetNorm = norm(targetY(:));
if targetNorm == 0
    targetNorm = 1;
end

bestShift = 0;
bestLoss = Inf;
for shiftIdx = 0:(nPoints - 1)
    shiftedY = circshift(candidateY, shiftIdx, 1);
    shiftLoss = norm(shiftedY(:) - targetY(:)) / targetNorm;
    if shiftLoss < bestLoss
        bestLoss = shiftLoss;
        bestShift = shiftIdx;
    end
end

alignedCandidateY = circshift(candidateY, bestShift, 1);
end
