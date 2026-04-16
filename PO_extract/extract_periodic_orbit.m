function result = extract_periodic_orbit(odefun, y0, parameter, opts)
%EXTRACT_PERIODIC_ORBIT Unified periodic-orbit solver entrypoint.
% Formal I/O contract: see PO_extract/IO_CONTRACT.md.

if nargin < 4 || isempty(opts)
    opts = struct();
end

result = po_initialize_result(y0, "direct");
if ~isstruct(opts)
    result.status = "invalid_options";
    result.message = "opts must be a struct.";
    return;
end

[opts, orbitSolver, message] = normalize_router_options(opts);
if strlength(message) > 0
    result.status = "invalid_options";
    result.message = message;
    result.diagnostics.requestedSolver = orbitSolver;
    return;
end

switch orbitSolver
    case "direct"
        result = solve_periodic_orbit_direct(odefun, y0, parameter, opts);
    case "matcont"
        result = run_matcont(odefun, y0, parameter, opts);
    otherwise
        result.status = "invalid_options";
        result.message = "Unsupported opts.orbit_solver value.";
end

result.diagnostics.requestedSolver = orbitSolver;
end

function result = run_matcont(odefun, y0, parameter, opts)
try
    result = solve_periodic_orbit_matcont(odefun, y0, parameter, opts);
catch ME
    result = po_initialize_result(y0, "matcont");
    result.status = "detection_failed";
    result.message = "MATCONT periodic-orbit solve failed: " + string(ME.message);
    result.diagnostics.matcont = struct( ...
        'errorIdentifier', string(ME.identifier), ...
        'errorMessage', string(ME.message));
    result = po_maybe_attach_quality(result, [], is_quality_requested(opts));
end
end

function [opts, orbitSolver, message] = normalize_router_options(opts)
message = "";
orbitSolver = "direct";

if isfield(opts, 'orbit_solver') && ~isempty(opts.orbit_solver)
    orbitSolver = string(opts.orbit_solver);
end

validSolvers = ["direct", "matcont"];
if ~any(orbitSolver == validSolvers)
    message = "opts.orbit_solver must be 'direct' or 'matcont'.";
    return;
end

opts.orbit_solver = orbitSolver;
end

function enabled = is_quality_requested(opts)
enabled = isstruct(opts) && isfield(opts, 'computeQualityDiagnostics') && ...
    ~isempty(opts.computeQualityDiagnostics) && logical(opts.computeQualityDiagnostics);
end
