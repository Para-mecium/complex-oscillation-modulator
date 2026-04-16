function result = extract_periodic_orbit(odefun, y0, parameter, opts)
%EXTRACT_PERIODIC_ORBIT Refine a periodic orbit via MATCONT from time integration.
%
% This interface is intentionally independent from PO_extract/extract_periodic_orbit.
% It preserves the four-argument call shape, but its option contract is
% MATCONT-specific and narrower than the direct backend.

if nargin < 4 || isempty(opts)
    opts = struct();
end

opts = normalize_options(opts);
[solverHandle, solverName, y0, opts, matcontOdefile, matcontParams, rhs] = ...
    validate_inputs(odefun, y0, parameter, opts);

result = initialize_result(y0);
result.backend_used = "matcont_orbit";
result = initialize_parameter_tracking(result, matcontParams, opts);
targetUserfunction = build_target_userfunction_odefile(matcontOdefile, opts.matcont_active_parameter, result.input_active_parameter_value);
result.diagnostics = struct( ...
    'solverName', string(solverName), ...
    'matcontOdefile', string(func2str(matcontOdefile)), ...
    'activeParameter', reshape(double(opts.matcont_active_parameter), 1, []), ...
    'matcontRoot', string(opts.matcont_root), ...
    'periodEstimate', NaN, ...
    'cycleWindowDuration', NaN, ...
    'tailStartTime', NaN, ...
    'refinedColumns', 0, ...
    'matcontStatus', "", ...
    'targetActiveParameter', result.input_active_parameter_value, ...
    'targetUserfunctionLabel', string(targetUserfunction.userInfo.label), ...
    'returnScanDirections', strings(0, 1), ...
    'returnScanColumns', zeros(0, 1), ...
    'returnScanHitBudget', false(0, 1));

ensure_matcont_paths(opts.matcont_root);

[tPilot, yPilot] = solverHandle(rhs, opts.tspan, y0, opts.odeOptions);

result.raw.pilot = struct('t', tPilot, 'y', yPilot);

[tTail, yTail, tTransient] = extract_tail_segment(tPilot, yPilot, opts);
result.diagnostics.tailStartTime = tTransient;

[periodEstimate, estimateDiagnostics] = estimate_period_from_tail(tTail, yTail);
result.diagnostics.periodEstimate = periodEstimate;
result.raw.period_estimate = estimateDiagnostics;
if ~isempty(opts.matcont_window_timespan)
    cycleWindowDuration = opts.matcont_window_timespan;
elseif isfinite(periodEstimate) && periodEstimate > 0
    cycleWindowDuration = opts.matcont_cycle_window_factor * periodEstimate;
else
    error('orbitmatcont:NoPositivePeriodEstimate', ...
        'Could not estimate a positive period from the pilot trajectory tail.');
end
result.diagnostics.cycleWindowDuration = cycleWindowDuration;

[tWindow, yWindow] = solverHandle(rhs, [0, cycleWindowDuration], yPilot(end, :).', opts.odeOptions);
result.raw.window = struct('t', tWindow, 'y', yWindow);

global cds lds
cds = [];
lds = [];

[x0, v0] = initOrbLC(targetUserfunction.odefile, tWindow, yWindow, matcontParams(:), ...
    opts.matcont_active_parameter, opts.matcont_ntst, opts.matcont_ncol, opts.matcont_tolerance);

if isempty(x0)
    error('orbitmatcont:InitOrbLCNoCycle', ...
        'MATCONT initOrbLC did not find a cycle in the selected orbit window.');
end

result.raw.seed = struct('x0', x0, 'v0', v0);

seedOrbit = build_orbit_from_lc_column(x0, lds, opts.extractNumPoints);
result.raw.seedOrbit = seedOrbit;

opt = build_continuation_options(opts.matcont_options, struct('MaxNumPoints', 1, 'Backward', 0));
[xLC, vLC, sLC, hLC, fLC] = cont(@limitcycle, x0, v0, opt); %#ok<ASGLU>
result.diagnostics.refinedColumns = size(xLC, 2);
result.diagnostics.matcontStatus = "continuation_success";

result.raw.matcont = struct('x', xLC, 'v', vLC, 's', sLC);
correctedColumn = xLC(:, 1);
correctedTangent = [];
if ~isempty(vLC)
    correctedTangent = vLC(:, 1);
end
result.seed_corrected_parameter_value = get_active_parameter_value(correctedColumn, lds);

returnInfo = return_to_input_parameter(correctedColumn, correctedTangent, opts, targetUserfunction.userInfo);
result.diagnostics.returnScanDirections = string(returnInfo.scanDirections(:));
result.diagnostics.returnScanColumns = reshape(double(returnInfo.scanColumns), [], 1);
result.diagnostics.returnScanHitBudget = reshape(logical(returnInfo.scanHitBudget), [], 1);
result.raw.parameter_return = returnInfo.raw;

result.output_parameter_values = reshape(extract_parameter_values(returnInfo.correctedColumn, lds), 1, []);
result.output_active_parameter_value = get_active_parameter_value(returnInfo.correctedColumn, lds);
result.parameter_error = abs(result.output_active_parameter_value - result.input_active_parameter_value);
result.parameter_status = string(returnInfo.parameterStatus);

orbit = build_orbit_from_lc_column(returnInfo.correctedColumn, lds, opts.extractNumPoints);
result.success = true;
result.has_orbit = true;
result.status = "converged_periodic_orbit";
result.code = 2;
result.message = "Periodic orbit initialized from time integration, corrected onto the LC branch, and matched back to the input parameter.";
result.period = orbit.period;
result.orbit_t = orbit.t;
result.orbit_y = orbit.y;
[result.max_variable, result.min_variable, result.amplitude] = compute_statistics(orbit.y);
end

function opts = normalize_options(opts)
cfg = orbitmatcont.default_config();
defaults = struct( ...
    'solver', [], ...
    'solver_name', cfg.solver_name, ...
    'odeOptions', struct(), ...
    'solver_options', struct(), ...
    'solver_tol', cfg.solver_tol, ...
    'tspan', [], ...
    'total_timespan', [], ...
    'single_timespan', [], ...
    'max_windows', [], ...
    'transientFraction', cfg.transientFraction, ...
    'transientTime', [], ...
    'extractNumPoints', cfg.extractNumPoints, ...
    'matcont_root', cfg.matcont_root, ...
    'matcont_odefile', cfg.matcont_odefile, ...
    'matcont_active_parameter', cfg.matcont_active_parameter, ...
    'matcont_parameter_values', cfg.matcont_parameter_values, ...
    'matcont_ntst', cfg.matcont_ntst, ...
    'matcont_ncol', cfg.matcont_ncol, ...
    'matcont_tolerance', cfg.matcont_tolerance, ...
    'matcont_window_timespan', cfg.matcont_window_timespan, ...
    'matcont_cycle_window_factor', cfg.matcont_cycle_window_factor, ...
    'matcont_options', cfg.matcont_options, ...
    'matcont_parameter_tolerance', cfg.matcont_parameter_tolerance, ...
    'matcont_return_max_points', cfg.matcont_return_max_points, ...
    'matcont_return_scan_both_directions', cfg.matcont_return_scan_both_directions);

names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
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

if isempty(opts.tspan)
    if ~isempty(opts.total_timespan)
        opts.tspan = [0, opts.total_timespan];
    elseif ~isempty(opts.single_timespan)
        maxWindows = max(1, double(default_scalar(opts.max_windows, 1)));
        opts.tspan = [0, opts.single_timespan * maxWindows];
    else
        error('orbitmatcont:MissingTspan', ...
            'Provide opts.tspan, opts.total_timespan, or opts.single_timespan.');
    end
end

baseOptions = opts.solver_options;
if isempty(baseOptions)
    baseOptions = struct();
end
if isstruct(baseOptions) && isstruct(opts.odeOptions) && ~isempty(opts.odeOptions)
    baseOptions = odeset(baseOptions, opts.odeOptions);
end
if isstruct(baseOptions)
    opts.odeOptions = odeset(baseOptions, ...
        'RelTol', opts.solver_tol.RelTol, ...
        'AbsTol', opts.solver_tol.AbsTol);
else
    error('orbitmatcont:InvalidSolverOptions', ...
        'opts.solver_options must be empty or a struct created by odeset.');
end
end

function [solverHandle, solverName, y0, opts, matcontOdefile, matcontParams, rhs] = ...
        validate_inputs(odefun, y0, parameter, opts)
solverHandle = [];
solverName = "";
matcontOdefile = [];
matcontParams = [];
rhs = [];
paramsForPilot = parameter;
matcontParams = extract_matcont_parameters(parameter, opts);

y0 = y0(:);
if isempty(y0) || ~isnumeric(y0) || ~all(isfinite(y0))
    error('orbitmatcont:InvalidInitialState', ...
        'Initial state y0 must be a finite numeric vector.');
end

if ~isstruct(opts.solver_tol) || ~isfield(opts.solver_tol, 'RelTol') || ~isfield(opts.solver_tol, 'AbsTol')
    error('orbitmatcont:InvalidSolverTolerance', ...
        'opts.solver_tol must contain RelTol and AbsTol.');
end

solverName = char(string(opts.solver_name));
if isa(opts.solver, 'function_handle')
    solverHandle = opts.solver;
elseif ischar(opts.solver) || isstring(opts.solver)
    solverHandle = str2func(char(string(opts.solver)));
else
    error('orbitmatcont:InvalidSolver', ...
        'opts.solver must be a solver name or function handle.');
end

if ~(isscalar(opts.matcont_active_parameter) && isnumeric(opts.matcont_active_parameter) && ...
        isfinite(opts.matcont_active_parameter) && opts.matcont_active_parameter >= 1 && ...
        abs(opts.matcont_active_parameter - round(opts.matcont_active_parameter)) < eps)
    error('orbitmatcont:InvalidActiveParameter', ...
        'opts.matcont_active_parameter must be a positive integer scalar.');
end

if ~(isscalar(opts.matcont_ntst) && isnumeric(opts.matcont_ntst) && isfinite(opts.matcont_ntst) && opts.matcont_ntst >= 2)
    error('orbitmatcont:InvalidNtst', ...
        'opts.matcont_ntst must be a finite integer >= 2.');
end

if ~(isscalar(opts.matcont_ncol) && isnumeric(opts.matcont_ncol) && isfinite(opts.matcont_ncol) && opts.matcont_ncol >= 2)
    error('orbitmatcont:InvalidNcol', ...
        'opts.matcont_ncol must be a finite integer >= 2.');
end

if abs(opts.matcont_ncol - round(opts.matcont_ncol)) >= eps || opts.matcont_ncol > 7
    error('orbitmatcont:InvalidNcol', ...
        'opts.matcont_ncol must be an integer in the MATCONT-supported range 2 <= ncol <= 7.');
end

if ~(isscalar(opts.matcont_tolerance) && isnumeric(opts.matcont_tolerance) && isfinite(opts.matcont_tolerance) && opts.matcont_tolerance > 0)
    error('orbitmatcont:InvalidMatcontTolerance', ...
        'opts.matcont_tolerance must be a positive finite scalar.');
end

if ~isempty(opts.matcont_window_timespan) && ...
        ~(isscalar(opts.matcont_window_timespan) && isnumeric(opts.matcont_window_timespan) && ...
        isfinite(opts.matcont_window_timespan) && opts.matcont_window_timespan > 0)
    error('orbitmatcont:InvalidWindowTimespan', ...
        'opts.matcont_window_timespan must be empty or a positive finite scalar.');
end

if ~(isscalar(opts.matcont_cycle_window_factor) && isnumeric(opts.matcont_cycle_window_factor) && ...
        isfinite(opts.matcont_cycle_window_factor) && opts.matcont_cycle_window_factor > 1 && ...
        opts.matcont_cycle_window_factor < 2)
    error('orbitmatcont:InvalidWindowFactor', ...
        'opts.matcont_cycle_window_factor must lie in (1, 2).');
end

if ~(isscalar(opts.matcont_parameter_tolerance) && isnumeric(opts.matcont_parameter_tolerance) && ...
        isfinite(opts.matcont_parameter_tolerance) && opts.matcont_parameter_tolerance > 0)
    error('orbitmatcont:InvalidParameterTolerance', ...
        'opts.matcont_parameter_tolerance must be a positive finite scalar.');
end

if ~(isscalar(opts.matcont_return_max_points) && isnumeric(opts.matcont_return_max_points) && ...
        isfinite(opts.matcont_return_max_points) && opts.matcont_return_max_points >= 2)
    error('orbitmatcont:InvalidReturnMaxPoints', ...
        'opts.matcont_return_max_points must be a finite integer >= 2.');
end

if ~(isscalar(opts.matcont_return_scan_both_directions) || isnumeric(opts.matcont_return_scan_both_directions))
    error('orbitmatcont:InvalidReturnDirectionsFlag', ...
        'opts.matcont_return_scan_both_directions must be a scalar logical.');
end

if isempty(matcontParams) || ~isnumeric(matcontParams) || ~isvector(matcontParams) || ~all(isfinite(matcontParams))
    error('orbitmatcont:InvalidParameterVector', ['MATCONT requires a finite numeric parameter vector. Provide numeric "parameter" ' ...
        'or opts.matcont_parameter_values.']);
end
matcontParams = matcontParams(:);
activeIndex = round(opts.matcont_active_parameter);
if activeIndex > numel(matcontParams)
    error('orbitmatcont:ActiveParameterOutOfRange', ...
        'opts.matcont_active_parameter exceeds the number of MATCONT parameter values.');
end

[matcontOdefile, rhs, paramsForPilot] = resolve_system_handles(odefun, paramsForPilot, opts, matcontParams);

end

function [matcontOdefile, rhs, parameterForPilot] = resolve_system_handles(odefun, parameterForPilot, opts, matcontParams)
matcontOdefile = [];
rhs = [];

if ~isempty(opts.matcont_odefile)
    matcontOdefile = opts.matcont_odefile;
elseif isa(odefun, 'function_handle') && nargin(odefun) == 0
    matcontOdefile = odefun;
else
    error('orbitmatcont:MissingMatcontOdefile', ['Provide opts.matcont_odefile, or pass a MATCONT odefile as the first argument ' ...
        'when using orbitmatcont.extract_periodic_orbit.']);
end

if ~isa(matcontOdefile, 'function_handle')
    error('orbitmatcont:InvalidMatcontOdefile', ...
        'opts.matcont_odefile must be a function handle.');
end

if isequal(matcontOdefile, odefun)
    parameterForPilot = matcontParams;
    handles = feval(matcontOdefile);
    if ~(iscell(handles) && numel(handles) >= 2)
        error('orbitmatcont:InvalidMatcontOdefile', ...
            'The MATCONT odefile must return at least the init and RHS handles.');
    end
    rhsFun = handles{2};
    rhs = @(t, y) call_matcont_rhs(rhsFun, t, y, matcontParams);
else
    if ~isa(odefun, 'function_handle')
        error('orbitmatcont:InvalidOdeFunction', ...
            'ODE function must be a function handle.');
    end
    rhs = @(t, y) call_odefun(odefun, t, y, parameterForPilot);
end
end

function params = extract_matcont_parameters(parameter, opts)
if isfield(opts, 'matcont_parameter_values') && ~isempty(opts.matcont_parameter_values)
    params = opts.matcont_parameter_values;
else
    params = parameter;
end
end

function ensure_matcont_paths(matcontRoot)
persistent configuredRoots
if isempty(configuredRoots)
    configuredRoots = strings(0, 1);
end

root = string(matcontRoot);
if any(configuredRoots == root)
    return;
end
if exist(matcontRoot, 'dir') ~= 7
    error('MATCONT root directory was not found: %s', matcontRoot);
end
addpath(genpath(matcontRoot));
configuredRoots(end + 1, 1) = root;
end

function [tTail, yTail, tTransient] = extract_tail_segment(t, y, opts)
t = t(:);
if isempty(opts.transientTime)
    tTransient = t(1) + opts.transientFraction * (t(end) - t(1));
else
    tTransient = opts.transientTime;
end
mask = t >= tTransient;
tTail = t(mask);
yTail = y(mask, :);
if numel(tTail) < 10
    error('Transient cutoff leaves too little tail data for MATCONT initialization.');
end
end

function [periodEstimate, diagnostics] = estimate_period_from_tail(tTail, yTail)
yTail = double(yTail);
anchor = yTail(1, :).';
searchStart = max(2, round(size(yTail, 1) / 3));
tailRange = max(yTail, [], 1) - min(yTail, [], 1);
tailScale = max(tailRange, 1e-12);
scaled = (yTail - anchor.') / tailScale;
distance = sum(abs(scaled(searchStart:end, :)), 2);
[bestDistance, idx] = min(distance);
bestIndex = searchStart + idx - 1;
periodEstimate = tTail(bestIndex) - tTail(1);
diagnostics = struct( ...
    'anchorTime', tTail(1), ...
    'searchStartTime', tTail(searchStart), ...
    'bestDistance', bestDistance, ...
    'bestIndex', bestIndex);
end

function opt = build_continuation_options(userOptions, overrides)
opt = contset;
opt = contset(opt, 'Singularities', 0);
opt = contset(opt, 'Adapt', 5);
opt = apply_contset_struct(opt, userOptions);
if nargin >= 2 && ~isempty(overrides)
    overrideNames = fieldnames(overrides);
    for i = 1:numel(overrideNames)
        opt = contset(opt, overrideNames{i}, overrides.(overrideNames{i}));
    end
end
end

function opt = build_target_continuation_options(userOptions, targetUserInfo, parameterTolerance, overrides)
opt = build_continuation_options(userOptions, overrides);
opt = contset(opt, 'Userfunctions', 1);
opt = contset(opt, 'UserfunctionsInfo', targetUserInfo);
currentTolerance = contget(opt, 'TestTolerance', parameterTolerance);
opt = contset(opt, 'TestTolerance', min(currentTolerance, parameterTolerance));
end

function targetUserfunction = build_target_userfunction_odefile(baseOdefile, activeIndex, targetValue)
targetUserfunction = struct();
targetUserfunction.odefile = @wrapped_odefile;
targetUserfunction.userInfo = struct( ...
    'name', 'active_parameter_target', ...
    'state', 1, ...
    'label', 'PA');

    function out = wrapped_odefile
        baseHandles = feval(baseOdefile);
        nBase = numel(baseHandles);
        nOut = max(10, nBase + 1);
        out = cell(1, nOut);
        out(1:min(9, nBase)) = baseHandles(1:min(9, nBase));
        out{10} = @target_active_parameter;
        if nBase >= 10
            out(11:nBase + 1) = baseHandles(10:end);
        end
    end

    function value = target_active_parameter(~, ~, varargin)
        params = reshape(cell2mat(varargin(:).'), [], 1);
        value = params(activeIndex) - targetValue;
    end
end

function opt = apply_contset_struct(opt, userOptions)
if isempty(userOptions)
    return;
end
if ~isstruct(userOptions)
    error('opts.matcont_options must be empty or a struct.');
end
names = fieldnames(userOptions);
for i = 1:numel(names)
    opt = contset(opt, names{i}, userOptions.(names{i}));
end
end

function orbit = build_orbit_from_lc_column(xColumn, lds, numPoints)
ncoords = lds.ncoords;
nphase = lds.nphase;
tps = lds.tps;
period = xColumn(lds.PeriodIdx);
rawCycle = reshape(xColumn(1:ncoords), nphase, tps).';
rawTime = lds.finemsh(:) * period;

if nargin < 3 || isempty(numPoints)
    numPoints = tps;
end
numPoints = max(2, round(numPoints));
targetTime = linspace(0, period, numPoints).';
targetCycle = interp1(rawTime, rawCycle, targetTime, 'pchip');

orbit = struct( ...
    'period', period, ...
    't', targetTime, ...
    'y', targetCycle, ...
    'raw_t', rawTime, ...
    'raw_y', rawCycle);
end

function returnInfo = return_to_input_parameter(correctedColumn, correctedTangent, opts, targetUserInfo)
global lds

targetValue = get_active_parameter_value_from_vector(lds.P0, opts.matcont_active_parameter);
correctedValue = get_active_parameter_value(correctedColumn, lds);
tolerance = opts.matcont_parameter_tolerance;
returnInfo = struct( ...
    'correctedColumn', [], ...
    'parameterStatus', "", ...
    'scanDirections', {cell(0, 1)}, ...
    'scanColumns', zeros(0, 1), ...
    'scanHitBudget', false(0, 1), ...
    'raw', struct('targetValue', targetValue, 'forward', [], 'backward', []));

if abs(correctedValue - targetValue) <= tolerance
    returnInfo.correctedColumn = correctedColumn;
    returnInfo.parameterStatus = "already_at_target";
    return;
end

directions = {'forward'};
if logical(opts.matcont_return_scan_both_directions)
    directions{end + 1} = 'backward'; %#ok<AGROW>
end

for i = 1:numel(directions)
    direction = directions{i};
    [xScan, vScan, sScan] = scan_branch(correctedColumn, correctedTangent, opts, direction, targetUserInfo);
    returnInfo.scanDirections{end + 1, 1} = direction; %#ok<AGROW>
    if isempty(xScan)
        error('orbitmatcont:EmptyTargetContinuation', ...
            'MATCONT returned no LC columns during target continuation.');
    end

    paramValues = extract_active_parameter_series(xScan, lds);
    returnInfo.scanColumns(end + 1, 1) = size(xScan, 2); %#ok<AGROW>
    hitBudget = scan_hit_max_points(xScan, opts.matcont_return_max_points);
    returnInfo.scanHitBudget(end + 1, 1) = hitBudget; %#ok<AGROW>
    returnInfo.raw.(direction) = struct( ...
        'status', "scan_success", ...
        'values', paramValues, ...
        'x', xScan, ...
        'v', vScan, ...
        's', sScan, ...
        'hitMaxPoints', hitBudget);
    if hitBudget
        report_target_scan_budget(direction, opts.matcont_return_max_points, paramValues, targetValue);
    end

    [locatedColumn, locatedIndex] = locate_userfunction_target(xScan, sScan, targetUserInfo, targetValue, lds, tolerance);
    if ~isempty(locatedColumn)
        returnInfo.correctedColumn = locatedColumn;
        returnInfo.raw.targetPoint = struct( ...
            'direction', string(direction), ...
            'index', locatedIndex, ...
            'parameterValue', get_active_parameter_value(locatedColumn, lds));
        returnInfo.parameterStatus = "returned_to_target";
        return;
    end
end

error('orbitmatcont:ParameterReturnFailed', ...
    'MATCONT continued the LC branch but did not locate the target active-parameter value via userfunction zero detection.');
end

function [xScan, vScan, sScan] = scan_branch(x0, v0, opts, direction, targetUserInfo)
xScan = [];
vScan = [];
sScan = [];
isBackward = strcmp(direction, 'backward');
opt = build_target_continuation_options(opts.matcont_options, targetUserInfo, opts.matcont_parameter_tolerance, struct( ...
    'MaxNumPoints', round(opts.matcont_return_max_points), ...
    'Backward', double(isBackward)));
[xScan, vScan, sScan] = cont(@limitcycle, x0, v0, opt);
end

function hitBudget = scan_hit_max_points(xColumns, maxPoints)
hitBudget = size(xColumns, 2) >= round(maxPoints);
end

function report_target_scan_budget(direction, maxPoints, paramValues, targetValue)
paramMin = min(paramValues);
paramMax = max(paramValues);
fprintf(['[orbitmatcont] Target continuation (%s) reached MaxNumPoints=%d ', ...
    'before exiting the scan; active-parameter range [%0.16g, %0.16g], target=%0.16g.\n'], ...
    direction, round(maxPoints), paramMin, paramMax, targetValue);
end

function values = extract_active_parameter_series(xColumns, lds)
paramIndex = get_active_parameter_column_index(lds);
values = reshape(xColumns(paramIndex, :), 1, []);
end

function [targetColumn, targetIndex] = locate_userfunction_target(xColumns, sPoints, targetUserInfo, targetValue, lds, tolerance)
targetColumn = [];
targetIndex = [];
if isempty(sPoints)
    return;
end

matchMask = false(size(sPoints));
for i = 1:numel(sPoints)
    label = "";
    msg = "";
    if isfield(sPoints(i), 'label')
        label = string(sPoints(i).label);
    end
    if isfield(sPoints(i), 'msg')
        msg = string(sPoints(i).msg);
    end
    matchMask(i) = (label == string(targetUserInfo.label)) | (msg == string(targetUserInfo.name));
end

matching = find(matchMask);
if isempty(matching)
    return;
end

candidateIndices = zeros(numel(matching), 1);
candidateErrors = inf(numel(matching), 1);
for i = 1:numel(matching)
    pointIndex = sPoints(matching(i)).index;
    if isnumeric(pointIndex) && isscalar(pointIndex) && pointIndex >= 1 && pointIndex <= size(xColumns, 2)
        candidateIndices(i) = pointIndex;
        candidateErrors(i) = abs(get_active_parameter_value(xColumns(:, pointIndex), lds) - targetValue);
    end
end

valid = candidateIndices > 0;
if ~any(valid)
    return;
end

[bestError, bestIdx] = min(candidateErrors(valid));
indices = candidateIndices(valid);
targetIndex = indices(bestIdx);
if bestError <= tolerance
    targetColumn = xColumns(:, targetIndex);
else
    targetIndex = [];
end
end

function paramValue = get_active_parameter_value(xColumn, lds)
paramIndex = get_active_parameter_column_index(lds);
paramValue = xColumn(paramIndex);
end

function params = extract_parameter_values(xColumn, lds)
params = lds.P0(:);
activeValue = get_active_parameter_value(xColumn, lds);
params(lds.ActiveParams(1)) = activeValue;
end

function columnIndex = get_active_parameter_column_index(lds)
columnIndex = get_period_index(lds) + 1;
end

function periodIndex = get_period_index(lds)
if isfield(lds, 'PeriodIdx') && ~isempty(lds.PeriodIdx)
    periodIndex = lds.PeriodIdx;
else
    periodIndex = lds.ncoords + 1;
end
end

function value = get_active_parameter_value_from_vector(parameterValues, activeIndex)
parameterValues = parameterValues(:);
value = parameterValues(round(activeIndex));
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

function dydt = call_matcont_rhs(rhsFun, t, y, parameterValues)
paramCell = num2cell(parameterValues(:).');
dydt = rhsFun(t, y, paramCell{:});
dydt = dydt(:);
end

function [maxVariable, minVariable, amplitude] = compute_statistics(orbitY)
if isempty(orbitY)
    maxVariable = [];
    minVariable = [];
    amplitude = [];
    return;
end
maxVariable = max(orbitY, [], 1);
minVariable = min(orbitY, [], 1);
amplitude = (maxVariable - minVariable) / 2;
end

function value = default_scalar(value, fallback)
if isempty(value)
    value = fallback;
end
end

function result = initialize_parameter_tracking(result, parameterValues, opts)
parameterValues = parameterValues(:);
activeIndex = round(opts.matcont_active_parameter);
result.input_parameter_values = reshape(parameterValues, 1, []);
result.output_parameter_values = reshape(parameterValues, 1, []);
result.active_parameter_index = activeIndex;
result.input_active_parameter_value = parameterValues(activeIndex);
result.seed_corrected_parameter_value = NaN;
result.output_active_parameter_value = NaN;
result.parameter_error = NaN;
result.parameter_status = "";
end

function result = initialize_result(y0)
nState = numel(y0);
result = struct( ...
    'success', false, ...
    'has_orbit', false, ...
    'status', "", ...
    'code', NaN, ...
    'message', "", ...
    'period', [], ...
    'orbit_t', [], ...
    'orbit_y', zeros(0, nState), ...
    'amplitude', [], ...
    'max_variable', [], ...
    'min_variable', [], ...
    'event_times', zeros(0, 1), ...
    'event_states', zeros(0, nState), ...
    'backend_used', "matcont_orbit", ...
    'input_parameter_values', zeros(1, 0), ...
    'output_parameter_values', zeros(1, 0), ...
    'active_parameter_index', [], ...
    'input_active_parameter_value', NaN, ...
    'seed_corrected_parameter_value', NaN, ...
    'output_active_parameter_value', NaN, ...
    'parameter_error', NaN, ...
    'parameter_status', "", ...
    'diagnostics', struct(), ...
    'raw', struct());
end
