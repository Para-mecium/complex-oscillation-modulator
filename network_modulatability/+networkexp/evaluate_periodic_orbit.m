function result = evaluate_periodic_orbit(odeFunc, y0, params, solverName, searchWindow, solverOptions)
%EVALUATE_PERIODIC_ORBIT Call extract_periodic_orbit via a normalized result contract.

if nargin < 6 || isempty(solverOptions)
    solverOptions = struct();
end
[solverOptions, detectionOptions] = networkexp.split_periodic_orbit_options(solverOptions);
solverOptions = normalize_solver_options(solverOptions);
solverOptions = filter_solver_options_for_solver(solverOptions, solverName);

packageDir = fileparts(mfilename('fullpath'));
networkDir = fileparts(packageDir);
fmamDir = fileparts(networkDir);
poExtractDir = fullfile(fmamDir, 'PO_extract');

result = struct( ...
    'success', false, ...
    'hasOrbit', false, ...
    'isCandidate', false, ...
    'status', 'missing_extractor', ...
    'orbitStatus', '', ...
    'orbitCode', NaN, ...
    'message', '', ...
    'TS', {{[], []}}, ...
    'observables', struct( ...
        'amplitude', [], ...
        'period', [], ...
        'maxVariable', [], ...
        'minVariable', []), ...
    'solverName', char(solverName), ...
    'searchWindow', searchWindow, ...
    'poExtractDir', poExtractDir, ...
    'errorIdentifier', '', ...
    'diagnostics', struct());

if exist('extract_periodic_orbit', 'file') ~= 2 && exist(poExtractDir, 'dir') == 7
    addpath(poExtractDir);
end

if exist('extract_periodic_orbit', 'file') ~= 2
    result.message = sprintf('extract_periodic_orbit.m was not found in %s.', poExtractDir);
    return
end

opts = struct( ...
    'solver_name', char(solverName), ...
    'tspan', [0, searchWindow], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-6), ...
    'solver_options', solverOptions);
opts = merge_struct_fields(opts, detectionOptions);

try
    extractResult = extract_periodic_orbit(odeFunc, y0, params, opts);
catch ME
    result.status = 'extract_error';
    result.message = ME.message;
    result.errorIdentifier = ME.identifier;
    return
end

result.orbitStatus = char(string(extractResult.status));
result.orbitCode = extractResult.code;
result.message = char(string(extractResult.message));
result.diagnostics = extractResult.diagnostics;
result.hasOrbit = isfield(extractResult, 'has_orbit') && extractResult.has_orbit;
result.isCandidate = (extractResult.code == 1);

if ~result.hasOrbit
    result.status = 'non_periodic';
    return
end

result.TS = {extractResult.orbit_t, extractResult.orbit_y};
result.observables = struct( ...
    'amplitude', reshape(extractResult.amplitude, 1, []), ...
    'period', extractResult.period, ...
    'maxVariable', reshape(extractResult.max_variable, 1, []), ...
    'minVariable', reshape(extractResult.min_variable, 1, []));

if extractResult.code == 2
    result.success = true;
    result.status = 'success';
elseif extractResult.code == 1
    result.status = 'candidate';
else
    result.status = 'non_periodic';
end
end

function target = merge_struct_fields(target, overrides)
if isempty(overrides)
    return
end

fieldNames = fieldnames(overrides);
for idxField = 1:numel(fieldNames)
    fieldName = fieldNames{idxField};
    target.(fieldName) = overrides.(fieldName);
end
end

function solverOptions = normalize_solver_options(solverOptions)
if isempty(solverOptions)
    solverOptions = struct();
elseif ~isstruct(solverOptions)
    error('networkexp:InvalidSolverOptions', ...
        'solverOptions must be empty or a struct created by odeset.');
end
end

function solverOptions = filter_solver_options_for_solver(solverOptions, solverName)
if is_stiff_solver_name(solverName)
    return
end

if isfield(solverOptions, 'Jacobian')
    solverOptions = rmfield(solverOptions, 'Jacobian');
end
if isfield(solverOptions, 'JPattern')
    solverOptions = rmfield(solverOptions, 'JPattern');
end
end

function tf = is_stiff_solver_name(solverName)
tf = ismember(lower(char(string(solverName))), {'ode15s', 'ode23s', 'ode23t', 'ode23tb'});
end
