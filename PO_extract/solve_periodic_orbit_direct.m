function result = solve_periodic_orbit_direct(odefun, y0, parameter, opts)
%SOLVE_PERIODIC_ORBIT_DIRECT Detect periodic-orbit candidates from forward integration.

if nargin < 4 || isempty(opts)
    opts = struct();
end

result = po_initialize_result(y0, "direct");
[opts, parse_message] = normalize_options(opts);
if strlength(string(parse_message)) > 0
    result.status = "invalid_options";
    result.message = string(parse_message);
    result = po_maybe_attach_quality(result, [], is_quality_enabled(opts));
    return;
end

[is_valid, validation_message, solver_handle, solver_name, y0, opts] = ...
    validate_inputs(odefun, y0, opts);
if ~is_valid
    result.status = "invalid_input";
    result.message = string(validation_message);
    result = po_maybe_attach_quality(result, [], opts.computeQualityDiagnostics);
    return;
end

try
    result = detect_periodic_orbit_core(odefun, solver_handle, solver_name, y0, parameter, opts);
catch ME
    result = po_initialize_result(y0, "direct");
    result.status = "detection_failed";
    result.message = "Periodic-orbit detection failed: " + string(ME.message);
    result.diagnostics.solverName = string(solver_name);
    result = po_maybe_attach_quality(result, [], opts.computeQualityDiagnostics);
end
end

function result = detect_periodic_orbit_core(odefun, solver_handle, solver_name, y0, parameter, opts)
rhs = @(t, y) call_odefun(odefun, t, y, parameter);
obsFcn = opts.observableFcn;
tspan = opts.tspan(:).';
t0 = tspan(1);
tf = tspan(end);

result = po_initialize_result(y0, "direct");
result.backend_used = "direct";
result.diagnostics.solverName = string(solver_name);
result.diagnostics.thresholds = build_detection_thresholds(opts);
result = po_maybe_attach_quality(result, [], opts.computeQualityDiagnostics);

try
    [t1, y1] = solver_handle(rhs, tspan, y0, opts.odeOptions);
catch solver_error
    result.status = "solver_failed";
    result.message = "ODE solver failed during pilot integration: " + string(solver_error.message);
    return;
end

obs1 = evalObservable(obsFcn, t1, y1, parameter, odefun);
if isempty(opts.transientTime)
    tTransient = t0 + opts.transientFraction * (tf - t0);
else
    tTransient = opts.transientTime;
end

tailMask = (t1 >= tTransient);
tTail = t1(tailMask);
yTail = y1(tailMask, :);
obsTail = obs1(tailMask);

sectionInfo = struct('mode', "", 'level', NaN, 'direction', opts.sectionDirection);
outRaw = struct();
outRaw.t = t1;
outRaw.y = y1;
outRaw.observable = obs1;
outRaw.tTransient = tTransient;
outRaw.tail = struct( ...
    't', tTail, ...
    'y', yTail, ...
    'observable', obsTail, ...
    'peaks', [], ...
    'peakTimes', [], ...
    'troughs', [], ...
    'troughTimes', []);

if numel(tTail) < 5
    result.status = "no_periodic_orbit_detected_on_tspan";
    result.code = 0;
    result.message = "Transient cutoff leaves too little data.";
    result.raw = outRaw;
    return;
end

[pks, pkLocs] = runFindPeaks(obsTail, tTail, opts);
[trsNeg, trLocs] = runFindPeaks(-obsTail, tTail, opts);
trs = -trsNeg;
tailRange = max(obsTail) - min(obsTail);
tailStd = std(obsTail);

if opts.autoSection
    c = median(obsTail);
    sectionFcn = @(t, y) call_observable(obsFcn, t, y, parameter, odefun) - c;
    sectionInfo.mode = "auto_observable_median";
    sectionInfo.level = c;
else
    if ~isempty(opts.sectionFcn)
        sectionFcn = @(t, y) call_section_function(opts.sectionFcn, t, y, parameter, odefun);
        sectionInfo.mode = "user_sectionFcn";
        sectionInfo.level = NaN;
    else
        c = opts.sectionLevel;
        sectionFcn = @(t, y) call_observable(obsFcn, t, y, parameter, odefun) - c;
        sectionInfo.mode = "user_observable_level";
        sectionInfo.level = c;
    end
end

eventOptions = odeset(opts.odeOptions, 'Events', @(t, y) section_event(t, y, sectionFcn, opts.sectionDirection));
try
    [t2, y2, te, ye] = solver_handle(rhs, tspan, y0, eventOptions);
catch solver_error
    result.status = "solver_failed";
    result.code = NaN;
    result.message = "ODE solver failed during event integration: " + string(solver_error.message);
    result.raw = outRaw;
    return;
end

obs2 = evalObservable(obsFcn, t2, y2, parameter, odefun);
keepEvent = (te >= tTransient);
te = te(keepEvent);
ye = ye(keepEvent, :);
[te, ye] = deduplicateEvents(te, ye, opts.minEventSeparation);

outRaw.t = t2;
outRaw.y = y2;
outRaw.observable = obs2;
outRaw.te = te;
outRaw.ye = ye;
outRaw.section = sectionInfo;
outRaw.tail.peaks = pks;
outRaw.tail.peakTimes = pkLocs;
outRaw.tail.troughs = trs;
outRaw.tail.troughTimes = trLocs;
result.raw = outRaw;

result.event_times = te;
result.event_states = ye;
result.diagnostics.tailRange = tailRange;
result.diagnostics.tailStd = tailStd;
result.diagnostics.numCrossings = numel(te);

if (tailRange <= opts.nonoscAmpTol && tailStd <= opts.nonoscStdTol) || ...
        ((numel(pks) < 2 || numel(trs) < 2) && tailRange <= 10 * opts.nonoscAmpTol)
    result.status = "decaying_to_equilibrium_or_nonoscillatory";
    result.code = -1;
    result.message = ['Tail oscillation is too small; trajectory looks more like ' ...
        'decay to an equilibrium or a nonoscillatory state.'];
    return;
end

if numel(te) < opts.minCrossings
    result.status = "no_periodic_orbit_detected_on_tspan";
    result.code = 0;
    result.message = ['Too few Poincare crossings after transient. ' ...
        'No convincing periodic-orbit evidence on the current time span.'];
    return;
end

nCycles = numel(te) - 1;
periods = diff(te);
amplitudes = nan(nCycles, 1);
stateCycleScales = nan(nCycles, 1);
for k = 1:nCycles
    [tCycle, yCycle] = extract_cycle_from_timeseries(t2, y2, te(k), te(k + 1), opts.samplesPerCycle);
    obsCycle = evalObservable(obsFcn, tCycle, yCycle, parameter, odefun);
    amplitudes(k) = max(obsCycle) - min(obsCycle);
    stateCycleScales(k) = max(norm(max(yCycle, [], 1) - min(yCycle, [], 1), 2), 1);
end

if size(ye, 1) >= 2
    dY = diff(ye, 1, 1);
    poincareResidual = rowNorms(dY) ./ stateCycleScales;
else
    poincareResidual = zeros(0, 1);
end

if numel(periods) >= 2
    periodRelChange = abs(diff(periods)) ./ max(abs(periods(1:end-1)), eps);
else
    periodRelChange = zeros(0, 1);
end

if numel(amplitudes) >= 2
    ampBase = max(abs(amplitudes(1:end-1)), opts.nonoscAmpTol);
    amplitudeRelChange = abs(diff(amplitudes)) ./ ampBase;
else
    amplitudeRelChange = zeros(0, 1);
end

result.diagnostics.periods = periods;
result.diagnostics.amplitudes = amplitudes;
result.diagnostics.stateCycleScales = stateCycleScales;
result.diagnostics.poincareResidual = poincareResidual;
result.diagnostics.periodRelChange = periodRelChange;
result.diagnostics.amplitudeRelChange = amplitudeRelChange;

m = opts.consecutiveCycles;
converged = false;
if nCycles >= m
    idx = (nCycles - m + 1):(nCycles - 1);
    if ~isempty(idx)
        converged = all(poincareResidual(idx) <= opts.poincareTol) && ...
            all(periodRelChange(idx) <= opts.periodTol) && ...
            all(amplitudeRelChange(idx) <= opts.amplitudeTol);
    end
end

candidate = false;
if nCycles >= 3
    mLoose = min(max(3, m), nCycles);
    idxLoose = (nCycles - mLoose + 1):(nCycles - 1);
    recentPeriods = periods(end - mLoose + 1:end);
    recentAmps = amplitudes(end - mLoose + 1:end);
    cvPeriod = std(recentPeriods) / max(mean(abs(recentPeriods)), eps);
    cvAmp = std(recentAmps) / max(mean(abs(recentAmps)), opts.nonoscAmpTol);

    looseStable = all(poincareResidual(idxLoose) <= opts.poincareTolLoose) && ...
        all(periodRelChange(idxLoose) <= opts.periodTolLoose) && ...
        all(amplitudeRelChange(idxLoose) <= opts.amplitudeTolLoose);

    candidateResidualTol = max(10 * opts.poincareTolLoose, sqrt(max(opts.poincareTolLoose, eps)));
    residualTrending = ~isempty(poincareResidual) && (poincareResidual(end) <= candidateResidualTol);

    candidate = residualTrending && (median(abs(recentAmps)) > opts.nonoscAmpTol) && ...
        (looseStable || (cvPeriod <= opts.periodTolLoose && cvAmp <= opts.amplitudeTolLoose));

    result.diagnostics.cvPeriod_recent = cvPeriod;
    result.diagnostics.cvAmplitude_recent = cvAmp;
    result.diagnostics.candidateResidualTol = candidateResidualTol;
end

decayingLike = false;
if numel(amplitudes) >= 4
    recent = amplitudes(end-3:end);
    decayingLike = all(diff(recent) <= 0) && (recent(end) < 0.5 * recent(1)) && ...
        (recent(end) < 10 * opts.nonoscAmpTol);
end

if converged
    result.status = "converged_periodic_orbit";
    result.code = 2;
    result.message = ['Poincare residual, period stability, and amplitude stability ' ...
        'all pass strict thresholds.'];
elseif decayingLike
    result.status = "decaying_to_equilibrium_or_nonoscillatory";
    result.code = -1;
    result.message = ['Oscillation amplitude keeps shrinking in the tail; ' ...
        'trajectory looks more like decay than convergence to a stable periodic orbit.'];
elseif candidate
    result.status = "candidate_periodic_orbit_not_converged";
    result.code = 1;
    result.message = ['There is periodic-orbit evidence, but strict convergence ' ...
        'criteria are not yet satisfied on the current time span.'];
else
    result.status = "no_periodic_orbit_detected_on_tspan";
    result.code = 0;
    result.message = ['Poincare crossings exist, but period/amplitude/return-point ' ...
        'stability are not convincing on the current time span.'];
end

if any(result.code == [1, 2])
    orbitStruct = build_orbit_struct(t2, y2, obsFcn, odefun, parameter, te(end - 1), te(end), ye(end, :), opts.extractNumPoints);
    result.has_orbit = true;
    result.orbit = orbitStruct;
    result.period = orbitStruct.period;
    result.orbit_t = orbitStruct.t;
    result.orbit_y = orbitStruct.y;
    [result.max_variable, result.min_variable, result.amplitude] = po_compute_statistics(orbitStruct.y);
end

result.success = (result.code == 2);
result = po_maybe_attach_quality(result, rhs, opts.computeQualityDiagnostics);
end

function orbit = build_orbit_struct(t, y, obsFcn, odefun, parameter, tStart, tEnd, sectionPoint, numPoints)
[tCycleAbs, yCycle] = extract_cycle_from_timeseries(t, y, tStart, tEnd, numPoints);
obsCycle = evalObservable(obsFcn, tCycleAbs, yCycle, parameter, odefun);

orbit = struct( ...
    't', tCycleAbs - tCycleAbs(1), ...
    'y', yCycle, ...
    'observable', obsCycle, ...
    'period', tEnd - tStart, ...
    'sectionPoint', reshape(sectionPoint, 1, []), ...
    'tAbs', tCycleAbs);
end

function [tq, yq] = extract_cycle_from_timeseries(t, y, ta, tb, numPoints)
[tUnique, yUnique] = deduplicateTimeseries(t, y);
tq = linspace(ta, tb, numPoints).';
yq = interp1(tUnique, yUnique, tq, 'pchip');
end

function [tUnique, yUnique] = deduplicateTimeseries(t, y)
[tUnique, keepIdx] = unique(t(:), 'stable');
yUnique = y(keepIdx, :);
end

function [value, isterminal, direction] = section_event(t, y, sectionFcn, sectionDirection)
value = sectionFcn(t, y);
isterminal = 0;
direction = sectionDirection;
end

function [opts, message] = normalize_options(opts)
message = "";

eventProvided = isfield(opts, 'event') && ~isempty(opts.event);
observableProvided = isfield(opts, 'observableFcn') && ~isempty(opts.observableFcn);
sectionProvided = (isfield(opts, 'sectionFcn') && ~isempty(opts.sectionFcn)) || ...
    (isfield(opts, 'sectionLevel') && ~isempty(opts.sectionLevel));
autoSectionProvided = isfield(opts, 'autoSection') && ~isempty(opts.autoSection);

defaults = struct( ...
    'solver', [], ...
    'solver_name', 'ode45', ...
    'odeOptions', struct(), ...
    'solver_options', struct(), ...
    'solver_tol', struct('RelTol', 1e-7, 'AbsTol', 1e-9), ...
    'tspan', [], ...
    'total_timespan', [], ...
    'single_timespan', [], ...
    'max_windows', [], ...
    'observableFcn', [], ...
    'autoSection', [], ...
    'sectionFcn', [], ...
    'sectionLevel', [], ...
    'sectionDirection', [], ...
    'transientFraction', 0.5, ...
    'transientTime', [], ...
    'minCrossings', [], ...
    'consecutiveCycles', 3, ...
    'poincareTol', 1e-2, ...
    'periodTol', 1e-2, ...
    'amplitudeTol', 1e-2, ...
    'poincareTolLoose', 1e-2, ...
    'periodTolLoose', 5e-2, ...
    'amplitudeTolLoose', 5e-2, ...
    'nonoscAmpTol', 1e-6, ...
    'nonoscStdTol', 1e-6, ...
    'minPeakProminence', [], ...
    'minPeakDistance', [], ...
    'samplesPerCycle', 300, ...
    'extractNumPoints', 400, ...
    'computeQualityDiagnostics', false, ...
    'minEventSeparation', [], ...
    'event', 1, ...
    'matcont_root', [], ...
    'matcont_odefile', [], ...
    'matcont_active_parameter', [], ...
    'matcont_ntst', [], ...
    'matcont_ncol', [], ...
    'matcont_tolerance', []);

unsupportedFields = po_find_unsupported_option_fields(opts);
if ~isempty(unsupportedFields)
    message = "Unsupported options in opts.";
    return;
end

option_names = fieldnames(defaults);
for idx = 1:numel(option_names)
    name = option_names{idx};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end

if isempty(opts.solver)
    opts.solver = opts.solver_name;
elseif isa(opts.solver, 'function_handle')
    opts.solver_name = func2str(opts.solver);
else
    opts.solver_name = char(string(opts.solver));
end

if isa(opts.solver_name, 'function_handle')
    opts.solver = opts.solver_name;
    opts.solver_name = func2str(opts.solver_name);
elseif isempty(opts.solver_name) && isa(opts.solver, 'function_handle')
    opts.solver_name = func2str(opts.solver);
end

if isempty(opts.tspan)
    if ~isempty(opts.total_timespan)
        totalTime = opts.total_timespan;
    elseif ~isempty(opts.single_timespan) && ~isempty(opts.max_windows)
        totalTime = opts.single_timespan * opts.max_windows;
    elseif ~isempty(opts.single_timespan)
        totalTime = opts.single_timespan;
    else
        totalTime = 100;
    end
    opts.tspan = [0, totalTime];
end

if isempty(opts.minCrossings)
    opts.minCrossings = 6;
end

if isempty(opts.observableFcn)
    opts.observableFcn = @(t, y) y(1);
end

if isempty(opts.autoSection)
    opts.autoSection = true;
end

if isempty(opts.sectionDirection)
    opts.sectionDirection = +1;
end

if eventProvided || ~(observableProvided || sectionProvided || autoSectionProvided)
    [opts, event_message] = apply_legacy_event_mapping(opts);
    if strlength(string(event_message)) > 0
        message = event_message;
        return;
    end
end

baseOptions = opts.solver_options;
if isempty(baseOptions)
    baseOptions = struct();
end
if isstruct(baseOptions) && isstruct(opts.odeOptions) && ~isempty(opts.odeOptions)
    baseOptions = odeset(baseOptions, opts.odeOptions);
end
if isstruct(baseOptions) && isstruct(opts.odeOptions)
    opts.odeOptions = odeset(baseOptions, ...
        'RelTol', opts.solver_tol.RelTol, ...
        'AbsTol', opts.solver_tol.AbsTol, ...
        'Events', []);
end

t0 = opts.tspan(1);
tf = opts.tspan(end);
if isempty(opts.minEventSeparation)
    opts.minEventSeparation = 1e-10 * max(tf - t0, 1);
end
end

function enabled = is_quality_enabled(opts)
enabled = isstruct(opts) && isfield(opts, 'computeQualityDiagnostics') && ...
    ~isempty(opts.computeQualityDiagnostics) && logical(opts.computeQualityDiagnostics);
end

function [opts, message] = apply_legacy_event_mapping(opts)
message = "";
eventOpt = opts.event;

if isempty(eventOpt)
    return;
end

if isnumeric(eventOpt) && isscalar(eventOpt)
    idx = double(eventOpt);
    opts.observableFcn = @(t, y) y(idx);
    opts.autoSection = true;
    if isempty(opts.sectionDirection)
        opts.sectionDirection = +1;
    end
    return;
end

if isa(eventOpt, 'function_handle')
    opts.sectionFcn = eventOpt;
    opts.autoSection = false;
    if isempty(opts.sectionDirection)
        opts.sectionDirection = 0;
    end
    return;
end

if ~isstruct(eventOpt)
    message = "opts.event must be empty, a state index, a struct, or a function handle.";
    return;
end

if isfield(eventOpt, 'direction') && ~isempty(eventOpt.direction)
    opts.sectionDirection = eventOpt.direction;
end

if isfield(eventOpt, 'handle') && isa(eventOpt.handle, 'function_handle')
    opts.sectionFcn = eventOpt.handle;
    opts.autoSection = false;
    if isempty(opts.sectionDirection)
        opts.sectionDirection = 0;
    end
    return;
end

if isfield(eventOpt, 'index')
    idx = eventOpt.index;
elseif isfield(eventOpt, 'state_index')
    idx = eventOpt.state_index;
else
    idx = 1;
end

if isfield(eventOpt, 'level') && ~isempty(eventOpt.level)
    opts.sectionLevel = eventOpt.level;
    opts.autoSection = false;
else
    opts.autoSection = true;
end
opts.observableFcn = @(t, y) y(idx);
end

function [is_valid, message, solver_handle, solver_name, y0, opts] = validate_inputs(odefun, y0, opts)
is_valid = false;
message = "";
solver_handle = [];
solver_name = "";
y0 = y0(:);

if ~isa(odefun, 'function_handle')
    message = "ODE function must be a function handle.";
    return;
end

if isempty(y0) || ~isnumeric(y0) || ~all(isfinite(y0))
    message = "Initial state y0 must be a finite numeric vector.";
    return;
end

if ~isstruct(opts.odeOptions) || ~isstruct(opts.solver_options)
    message = "opts.odeOptions and opts.solver_options must be structs created by odeset.";
    return;
end

if ~isstruct(opts.solver_tol) || ~isfield(opts.solver_tol, 'RelTol') || ~isfield(opts.solver_tol, 'AbsTol')
    message = "opts.solver_tol must contain RelTol and AbsTol.";
    return;
end

if ~(isscalar(opts.solver_tol.RelTol) && isnumeric(opts.solver_tol.RelTol) && ...
        isfinite(opts.solver_tol.RelTol) && opts.solver_tol.RelTol > 0)
    message = "opts.solver_tol.RelTol must be a positive finite scalar.";
    return;
end

if ~(isnumeric(opts.solver_tol.AbsTol) && all(isfinite(opts.solver_tol.AbsTol(:))) && ...
        all(opts.solver_tol.AbsTol(:) > 0))
    message = "opts.solver_tol.AbsTol must contain positive finite values.";
    return;
end

tspan = opts.tspan;
if ~(isnumeric(tspan) && isvector(tspan) && numel(tspan) >= 2 && all(isfinite(tspan)) && all(diff(tspan) > 0))
    message = "opts.tspan must be a strictly increasing numeric vector with at least two entries.";
    return;
end

if ~(isscalar(opts.minCrossings) && opts.minCrossings >= 2 && opts.minCrossings == floor(opts.minCrossings))
    message = "opts.minCrossings must be an integer greater than or equal to 2.";
    return;
end

if ~(isscalar(opts.consecutiveCycles) && opts.consecutiveCycles >= 2 && opts.consecutiveCycles == floor(opts.consecutiveCycles))
    message = "opts.consecutiveCycles must be an integer greater than or equal to 2.";
    return;
end

if ~(isscalar(opts.samplesPerCycle) && opts.samplesPerCycle >= 3 && opts.samplesPerCycle == floor(opts.samplesPerCycle))
    message = "opts.samplesPerCycle must be an integer greater than or equal to 3.";
    return;
end

if ~(isscalar(opts.extractNumPoints) && opts.extractNumPoints >= 3 && opts.extractNumPoints == floor(opts.extractNumPoints))
    message = "opts.extractNumPoints must be an integer greater than or equal to 3.";
    return;
end

tolNames = {'poincareTol', 'periodTol', 'amplitudeTol', ...
    'poincareTolLoose', 'periodTolLoose', 'amplitudeTolLoose', ...
    'nonoscAmpTol', 'nonoscStdTol', 'minEventSeparation'};
for k = 1:numel(tolNames)
    name = tolNames{k};
    value = opts.(name);
    if ~(isscalar(value) && isnumeric(value) && isfinite(value) && value >= 0)
        message = "opts." + string(name) + " must be a nonnegative finite scalar.";
        return;
    end
end

if ~(isscalar(opts.transientFraction) && isnumeric(opts.transientFraction) && ...
        isfinite(opts.transientFraction) && opts.transientFraction >= 0 && opts.transientFraction <= 1)
    message = "opts.transientFraction must be a finite scalar in [0, 1].";
    return;
end

if ~isempty(opts.transientTime)
    if ~(isscalar(opts.transientTime) && isnumeric(opts.transientTime) && isfinite(opts.transientTime) && ...
            opts.transientTime >= tspan(1) && opts.transientTime <= tspan(end))
        message = "opts.transientTime must lie inside opts.tspan.";
        return;
    end
end

if ~isa(opts.observableFcn, 'function_handle')
    message = "opts.observableFcn must be a function handle.";
    return;
end

if ~isempty(opts.sectionFcn) && ~isa(opts.sectionFcn, 'function_handle')
    message = "opts.sectionFcn must be empty or a function handle.";
    return;
end

if ~(isscalar(opts.sectionDirection) && isnumeric(opts.sectionDirection) && ...
        any(opts.sectionDirection == [-1, 0, 1]))
    message = "opts.sectionDirection must be -1, 0, or +1.";
    return;
end

if ~opts.autoSection && isempty(opts.sectionFcn) && isempty(opts.sectionLevel)
    message = "When autoSection is false, provide opts.sectionFcn or opts.sectionLevel.";
    return;
end

if ~isempty(opts.sectionLevel) && ~(isscalar(opts.sectionLevel) && isnumeric(opts.sectionLevel) && isfinite(opts.sectionLevel))
    message = "opts.sectionLevel must be a finite scalar when provided.";
    return;
end

if isempty(opts.solver)
    solver_name = char(string(opts.solver_name));
    if ~(exist(solver_name, 'file') == 2 || exist(solver_name, 'builtin') == 5)
        message = "Requested ODE solver '" + string(solver_name) + "' is not available.";
        return;
    end
    solver_handle = str2func(solver_name);
else
    if isa(opts.solver, 'function_handle')
        solver_handle = opts.solver;
        solver_name = func2str(opts.solver);
    else
        solver_name = char(string(opts.solver));
        if ~(exist(solver_name, 'file') == 2 || exist(solver_name, 'builtin') == 5)
            message = "Requested ODE solver '" + string(solver_name) + "' is not available.";
            return;
        end
        solver_handle = str2func(solver_name);
    end
end

opts.solver = solver_handle;
opts.solver_name = solver_name;
is_valid = true;
end

function obs = evalObservable(obsFcn, t, y, parameter, odefun)
n = numel(t);
obs = zeros(n, 1);
for i = 1:n
    obs(i) = call_observable(obsFcn, t(i), y(i, :).', parameter, odefun);
end
end

function value = call_observable(obsFcn, t, y, parameter, odefun)
value = call_user_function(obsFcn, t, y, parameter, odefun);
if ~(isscalar(value) && isnumeric(value) && isfinite(value))
    error('Observable function must return a finite scalar.');
end
end

function value = call_section_function(sectionFcn, t, y, parameter, odefun)
value = call_user_function(sectionFcn, t, y, parameter, odefun);
if ~(isscalar(value) && isnumeric(value) && isfinite(value))
    error('Section function must return a finite scalar.');
end
end

function out = call_user_function(handle, t, y, parameter, odefun)
narginHandle = nargin(handle);
if narginHandle == 2
    out = handle(t, y);
elseif narginHandle == 3
    out = handle(t, y, parameter);
else
    out = handle(t, y, parameter, odefun);
end
end

function dydt = call_odefun(odefun, t, y, parameter)
narginHandle = nargin(odefun);
if narginHandle == 2
    dydt = odefun(t, y);
else
    dydt = odefun(t, y, parameter);
end
dydt = dydt(:);
end

function [pks, locs] = runFindPeaks(sig, t, opts)
args = {};
if ~isempty(opts.minPeakProminence)
    args = [args, {'MinPeakProminence', opts.minPeakProminence}];
end
if ~isempty(opts.minPeakDistance)
    args = [args, {'MinPeakDistance', opts.minPeakDistance}];
end

try
    [pks, locs] = findpeaks(sig, t, args{:});
catch
    [pks, locs] = fallbackFindPeaks(sig, t);
end
end

function [pks, locs] = fallbackFindPeaks(sig, t)
sig = sig(:);
t = t(:);
if numel(sig) < 3
    pks = zeros(0, 1);
    locs = zeros(0, 1);
    return;
end

mask = (sig(2:end-1) >= sig(1:end-2)) & (sig(2:end-1) >= sig(3:end));
idx = find(mask) + 1;
pks = sig(idx);
locs = t(idx);
end

function [te2, ye2] = deduplicateEvents(te, ye, minSep)
if isempty(te)
    te2 = te;
    ye2 = ye;
    return;
end

keep = true(size(te));
for i = 2:numel(te)
    if abs(te(i) - te(i - 1)) <= minSep
        keep(i) = false;
    end
end
te2 = te(keep);
ye2 = ye(keep, :);
end

function r = rowNorms(A)
r = sqrt(sum(A.^2, 2));
end

function thresholds = build_detection_thresholds(opts)
thresholds = struct( ...
    'minCrossings', opts.minCrossings, ...
    'consecutiveCycles', opts.consecutiveCycles, ...
    'nonoscAmpTol', opts.nonoscAmpTol, ...
    'nonoscStdTol', opts.nonoscStdTol, ...
    'strict', struct( ...
        'poincareTol', opts.poincareTol, ...
        'periodTol', opts.periodTol, ...
        'amplitudeTol', opts.amplitudeTol), ...
    'loose', struct( ...
        'poincareTol', opts.poincareTolLoose, ...
        'periodTol', opts.periodTolLoose, ...
        'amplitudeTol', opts.amplitudeTolLoose));
end
