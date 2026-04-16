function features = evaluate_orbit_features(orbit, observableSpec, featureSpec, opts)
%EVALUATE_ORBIT_FEATURES Evaluate orbit-derived state and observable features.

if nargin < 2
    observableSpec = [];
end
if nargin < 3
    featureSpec = [];
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

opts = normalize_options(opts);
orbit = normalize_orbit_input(orbit, opts);
observableFns = normalize_observable_spec(observableSpec);
requestedFeatures = normalize_feature_spec(featureSpec, ~isempty(observableFns));

features = struct( ...
    'period', orbit.period, ...
    'state', struct('series', orbit.y, 'max', [], 'min', [], 'amplitude', []), ...
    'observable', struct('series', zeros(size(orbit.y, 1), 0), 'max', [], 'min', [], 'amplitude', []), ...
    'diagnostics', struct( ...
        'requested_features', requestedFeatures, ...
        'state', empty_group_diagnostics(0), ...
        'observable', empty_group_diagnostics(0)));

stateRequested = any(startsWith(requestedFeatures, "state_"));
observableRequested = any(startsWith(requestedFeatures, "observable_"));

if stateRequested
    [features.state, features.diagnostics.state] = compute_feature_group( ...
        orbit.t, orbit.y, orbit.period, opts);
end

if observableRequested
    if isempty(observableFns)
        error('orbitfeature:ObservableRequired', ...
            'observableSpec is required when requesting observable features.');
    end
    observableSeries = evaluate_observables(observableFns, orbit.t, orbit.y);
    [features.observable, features.diagnostics.observable] = compute_feature_group( ...
        orbit.t, observableSeries, orbit.period, opts);
end
end

function opts = normalize_options(opts)
if ~isstruct(opts)
    error('orbitfeature:InvalidOptions', 'opts must be a struct.');
end

supportedFields = ["extremaRefinement", "refinementFactor", "refinementPointCount", ...
    "refinementMethod", "periodTolerance", "endpointTolerance"];
optionNames = string(fieldnames(opts));
unsupported = sort(optionNames(~ismember(optionNames, supportedFields)));
if ~isempty(unsupported)
    error('orbitfeature:UnsupportedOption', ...
        'Unsupported options: %s.', strjoin(cellstr(unsupported), ', '));
end

defaults = struct( ...
    'extremaRefinement', true, ...
    'refinementFactor', 100, ...
    'refinementPointCount', 7, ...
    'refinementMethod', "spline", ...
    'periodTolerance', 1e-8, ...
    'endpointTolerance', 1e-8);

names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end

if ~(isscalar(opts.extremaRefinement) && ...
        (islogical(opts.extremaRefinement) || isnumeric(opts.extremaRefinement)))
    error('orbitfeature:InvalidOptions', ...
        'opts.extremaRefinement must be a logical scalar.');
end
opts.extremaRefinement = logical(opts.extremaRefinement);

validate_positive_integer_option(opts.refinementFactor, 'refinementFactor');
validate_odd_point_count_option(opts.refinementPointCount, 'refinementPointCount');
opts.refinementMethod = normalize_refinement_method(opts.refinementMethod);
validate_positive_scalar_option(opts.periodTolerance, 'periodTolerance');
validate_positive_scalar_option(opts.endpointTolerance, 'endpointTolerance');
end

function validate_positive_scalar_option(value, name)
if ~(isscalar(value) && isnumeric(value) && isfinite(value) && value > 0)
    error('orbitfeature:InvalidOptions', ...
        'opts.%s must be a positive finite scalar.', name);
end
end

function validate_positive_integer_option(value, name)
if ~(isscalar(value) && isnumeric(value) && isreal(value) && isfinite(value) && ...
        value >= 2 && value == floor(value))
    error('orbitfeature:InvalidOptions', ...
        'opts.%s must be an integer scalar greater than or equal to 2.', name);
end
end

function validate_odd_point_count_option(value, name)
if ~(isscalar(value) && isnumeric(value) && isreal(value) && isfinite(value) && ...
        value >= 3 && value == floor(value) && mod(value, 2) == 1)
    error('orbitfeature:InvalidOptions', ...
        'opts.%s must be an odd integer scalar greater than or equal to 3.', name);
end
end

function method = normalize_refinement_method(value)
if ischar(value) || (isstring(value) && isscalar(value))
    method = lower(string(value));
else
    error('orbitfeature:InvalidOptions', ...
        'opts.refinementMethod must be a string scalar or character vector.');
end

supportedMethods = ["quadratic", "spline"];
if ~ismember(method, supportedMethods)
    error('orbitfeature:InvalidOptions', ...
        'opts.refinementMethod must be one of: %s.', ...
        strjoin(cellstr(supportedMethods), ', '));
end
end

function orbit = normalize_orbit_input(orbit, opts)
requiredFields = ["t", "y", "period"];
if ~isstruct(orbit) || numel(orbit) ~= 1
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit must be a scalar struct with fields t, y, and period.');
end

orbitFields = string(fieldnames(orbit));
missing = requiredFields(~ismember(requiredFields, orbitFields));
if ~isempty(missing)
    error('orbitfeature:MissingOrbitField', ...
        'orbit is missing required fields: %s.', strjoin(cellstr(missing), ', '));
end

unsupported = sort(orbitFields(~ismember(orbitFields, requiredFields)));
if ~isempty(unsupported)
    error('orbitfeature:UnsupportedOrbitField', ...
        'Unsupported orbit fields: %s.', strjoin(cellstr(unsupported), ', '));
end

t = orbit.t(:);
y = orbit.y;
period = orbit.period;

if ~(isnumeric(t) && isreal(t) && ~isempty(t) && all(isfinite(t)))
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit.t must be a finite real numeric vector.');
end
if ~(isnumeric(y) && isreal(y) && ismatrix(y) && ~isempty(y) && all(isfinite(y(:))))
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit.y must be a finite real numeric matrix.');
end
if size(y, 1) ~= numel(t)
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit.y must have the same number of rows as orbit.t samples.');
end
if ~(isscalar(period) && isnumeric(period) && isreal(period) && isfinite(period) && period > 0)
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit.period must be a positive finite real scalar.');
end

t = t - t(1);
if any(diff(t) <= 0)
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit.t must be strictly increasing.');
end

periodTol = opts.periodTolerance * max(abs(period), 1);
if abs(t(end) - period) > periodTol
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit.t must end at orbit.period within tolerance.');
end

endpointScale = max(1, max(abs(y(:))));
closureTolerance = opts.endpointTolerance * endpointScale;
if size(y, 1) >= 2 && norm(y(end, :) - y(1, :), inf) <= closureTolerance
    t = t(1:end-1);
    y = y(1:end-1, :);
end

if isempty(t) || size(y, 1) < 2
    error('orbitfeature:InvalidOrbitInput', ...
        'orbit must contain at least two non-redundant samples over one period.');
end

orbit = struct( ...
    't', t, ...
    'y', y, ...
    'period', period);
end

function observableFns = normalize_observable_spec(observableSpec)
if isempty(observableSpec)
    observableFns = {};
    return;
end

if isa(observableSpec, 'function_handle')
    observableFns = {observableSpec};
    return;
end

if iscell(observableSpec)
    observableFns = observableSpec(:).';
    for i = 1:numel(observableFns)
        if ~isa(observableFns{i}, 'function_handle')
            error('orbitfeature:InvalidObservableSpec', ...
                'observableSpec cell entries must be function handles.');
        end
    end
    return;
end

error('orbitfeature:InvalidObservableSpec', ...
    'observableSpec must be empty, a function handle, or a cell array of function handles.');
end

function requestedFeatures = normalize_feature_spec(featureSpec, hasObservable)
supported = [ ...
    "period"
    "state_max"
    "state_min"
    "state_amplitude"
    "observable_max"
    "observable_min"
    "observable_amplitude"
    ];

if nargin < 2
    hasObservable = false;
end

if isempty(featureSpec)
    requestedFeatures = ["period", "state_max", "state_min", "state_amplitude"];
    if hasObservable
        requestedFeatures = [requestedFeatures, ...
            "observable_max", "observable_min", "observable_amplitude"];
    end
    return;
end

if ischar(featureSpec) || (isstring(featureSpec) && isscalar(featureSpec))
    requestedFeatures = string(featureSpec);
elseif iscellstr(featureSpec) || isstring(featureSpec) || iscell(featureSpec)
    requestedFeatures = string(featureSpec(:)).';
else
    error('orbitfeature:InvalidFeatureSpec', ...
        'featureSpec must be empty, a string, or a cell array of strings.');
end

unsupported = requestedFeatures(~ismember(requestedFeatures, supported));
if ~isempty(unsupported)
    error('orbitfeature:UnsupportedFeature', ...
        'Unsupported feature names: %s.', strjoin(cellstr(unique(unsupported, 'stable')), ', '));
end

requestedFeatures = unique(requestedFeatures, 'stable');
end

function [group, diagnostics] = compute_feature_group(t, series, period, opts)
series = double(series);
[maxValues, maxDiagnostics] = compute_extrema(t, series, period, "max", opts);
[minValues, minDiagnostics] = compute_extrema(t, series, period, "min", opts);

group = struct( ...
    'series', series, ...
    'max', maxValues, ...
    'min', minValues, ...
    'amplitude', 0.5 * (maxValues - minValues));

diagnostics = struct( ...
    'maxRefined', reshape([maxDiagnostics.refined], 1, []), ...
    'minRefined', reshape([minDiagnostics.refined], 1, []), ...
    'maxFallback', reshape([maxDiagnostics.fallback], 1, []), ...
    'minFallback', reshape([minDiagnostics.fallback], 1, []), ...
    'maxReason', string({maxDiagnostics.reason}), ...
    'minReason', string({minDiagnostics.reason}), ...
    'maxLocation', reshape([maxDiagnostics.location], 1, []), ...
    'minLocation', reshape([minDiagnostics.location], 1, []));
end

function diagnostics = empty_group_diagnostics(signalCount)
diagnostics = struct( ...
    'maxRefined', false(1, signalCount), ...
    'minRefined', false(1, signalCount), ...
    'maxFallback', false(1, signalCount), ...
    'minFallback', false(1, signalCount), ...
    'maxReason', strings(1, signalCount), ...
    'minReason', strings(1, signalCount), ...
    'maxLocation', nan(1, signalCount), ...
    'minLocation', nan(1, signalCount));
end

function [values, diagnostics] = compute_extrema(t, series, period, mode, opts)
nSignals = size(series, 2);
values = zeros(1, nSignals);
diagnostics = repmat(struct('refined', false, 'fallback', false, ...
    'reason', "", 'location', NaN), 1, nSignals);

for i = 1:nSignals
    [values(i), diagnostics(i)] = compute_single_extremum(t, series(:, i), period, mode, opts);
end
end

function [value, diagnostics] = compute_single_extremum(t, signal, period, mode, opts)
signal = signal(:);
diagnostics = struct('refined', false, 'fallback', false, 'reason', "", 'location', NaN);

switch mode
    case "max"
        [discreteValue, idx] = max(signal);
    case "min"
        [discreteValue, idx] = min(signal);
    otherwise
        error('orbitfeature:InternalError', 'Unsupported extremum mode.');
end

diagnostics.location = t(idx);

if ~opts.extremaRefinement
    diagnostics.fallback = true;
    diagnostics.reason = "refinement_disabled";
    value = discreteValue;
    return;
end

if numel(signal) < 3
    diagnostics.fallback = true;
    diagnostics.reason = "insufficient_samples";
    value = discreteValue;
    return;
end

if opts.refinementPointCount > numel(signal)
    diagnostics.fallback = true;
    diagnostics.reason = "insufficient_local_points";
    value = discreteValue;
    return;
end

try
    [xStar, refinedValue] = refine_extremum_on_dense_grid(t, signal, period, idx, ...
        mode, opts.refinementFactor, opts.refinementMethod, opts.refinementPointCount);
    if ~isfinite(refinedValue)
        error('orbitfeature:InvalidRefinementValue', ...
            'Refined extremum must be finite.');
    end
    value = refinedValue;
    diagnostics.refined = true;
    diagnostics.reason = "refined_dense_grid";
    diagnostics.location = xStar;
catch
    diagnostics.fallback = true;
    diagnostics.reason = "refinement_failed";
    value = discreteValue;
end
end

function [xStar, value] = refine_extremum_on_dense_grid(t, signal, period, idx, ...
    mode, refinementFactor, refinementMethod, refinementPointCount)
[denseT, denseSignal] = evaluate_refined_signal(t, signal, period, idx, refinementFactor, ...
    refinementMethod, refinementPointCount);
if any(~isfinite(denseSignal))
    error('orbitfeature:InvalidRefinementValue', ...
        'Dense-grid refinement values must be finite.');
end

switch mode
    case "max"
        [value, idx] = max(denseSignal);
    case "min"
        [value, idx] = min(denseSignal);
    otherwise
        error('orbitfeature:InternalError', 'Unsupported extremum mode.');
end

xStar = wrap_to_period(denseT(idx), period);
end

function [denseT, denseSignal] = evaluate_refined_signal(t, signal, period, idx, refinementFactor, ...
    method, refinementPointCount)
[localT, localSignal] = build_local_stencil(t, signal, period, idx, refinementPointCount);
switch method
    case "spline"
        denseCount = 2 * refinementFactor + 1;
        denseT = linspace(localT(1), localT(end), denseCount).';
        denseSignal = interp1(localT, localSignal, denseT, 'spline');
    case "quadratic"
        [denseT, denseSignal] = evaluate_quadratic_refinement(localT, localSignal, refinementFactor);
    otherwise
        error('orbitfeature:InternalError', 'Unsupported refinement method.');
end
end

function [denseT, denseSignal] = evaluate_quadratic_refinement(localT, localSignal, refinementFactor)
denseCount = 2 * refinementFactor + 1;
denseT = linspace(localT(1), localT(end), denseCount).';
centerT = localT((numel(localT) + 1) / 2);
coeffs = polyfit(localT - centerT, localSignal, 2);
denseSignal = polyval(coeffs, denseT - centerT);
end

function [localT, localSignal] = build_local_stencil(t, signal, period, idx, pointCount)
n = numel(t);
halfSpan = (pointCount - 1) / 2;
localT = zeros(pointCount, 1);
localSignal = zeros(pointCount, 1);
for offset = -halfSpan:halfSpan
    writeIdx = offset + halfSpan + 1;
    sampleIdx = idx + offset;
    sampleT = t(mod(sampleIdx - 1, n) + 1);
    if sampleIdx < 1
        sampleT = sampleT - period;
    elseif sampleIdx > n
        sampleT = sampleT + period;
    end
    localT(writeIdx) = sampleT;
    localSignal(writeIdx) = signal(mod(sampleIdx - 1, n) + 1);
end
end

function wrapped = wrap_to_period(x, period)
wrapped = mod(x, period);
if wrapped < 0
    wrapped = wrapped + period;
end
end

function observableSeries = evaluate_observables(observableFns, t, y)
observableSeries = zeros(numel(t), 0);
for i = 1:numel(observableFns)
    value = observableFns{i}(t, y);
    value = normalize_observable_output(value, numel(t), i);
    observableSeries = [observableSeries, value]; %#ok<AGROW>
end
end

function value = normalize_observable_output(value, nRows, idx)
if ~(isnumeric(value) && isreal(value) && all(isfinite(value(:))))
    error('orbitfeature:InvalidObservableOutput', ...
        'Observable %d must return a finite real numeric array.', idx);
end

if isvector(value) && numel(value) == nRows
    value = reshape(value, nRows, 1);
elseif size(value, 1) ~= nRows
    error('orbitfeature:InvalidObservableOutput', ...
        'Observable %d must return values with one row per time sample.', idx);
end
end
