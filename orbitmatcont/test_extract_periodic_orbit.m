function test_extract_periodic_orbit()
%TEST_EXTRACT_PERIODIC_ORBIT Regression checks for orbitmatcont.extract_periodic_orbit.

rootDir = fileparts(fileparts(mfilename('fullpath')));
matcontRoot = fullfile(rootDir, 'MatCont7p6');
addpath(rootDir);
addpath(genpath(matcontRoot));

opts = struct( ...
    'tspan', [0, 300], ...
    'transientFraction', 0.8, ...
    'matcont_root', matcontRoot, ...
    'matcont_active_parameter', 2, ...
    'matcont_tolerance', 1e-2, ...
    'matcont_parameter_tolerance', 1e-3, ...
    'matcont_window_timespan', 10, ...
    'extractNumPoints', 120);

directResult = orbitmatcont.extract_periodic_orbit(@adaptx, [0.3; 0.5; -0.1], [1; 0.8], opts);
assert(directResult.has_orbit, 'Expected orbit extraction to succeed for MATCONT odefile input.');
assert(strcmp(char(directResult.status), 'converged_periodic_orbit'), ...
    'Expected converged_periodic_orbit status.');
assert(directResult.success, 'Expected successful extraction for MATCONT odefile input.');
assert(directResult.period > 5 && directResult.period < 9, ...
    'Expected refined period to stay in a reasonable range.');
assert(size(directResult.orbit_y, 1) == opts.extractNumPoints, ...
    'Expected extracted orbit to be resampled to extractNumPoints.');
assert(abs(directResult.output_active_parameter_value - directResult.input_active_parameter_value) <= opts.matcont_parameter_tolerance, ...
    'Expected the returned active parameter to match the input parameter within tolerance.');
assert(~isempty(directResult.output_parameter_values), ...
    'Expected output_parameter_values to be populated.');
assert(strlength(string(directResult.parameter_status)) > 0, ...
    'Expected parameter_status to be reported.');

handles = adaptx();
rhs = @(t, y, p) handles{2}(t, y, p(1), p(2));
wrappedOpts = opts;
wrappedOpts.matcont_odefile = @adaptx;
wrappedResult = orbitmatcont.extract_periodic_orbit(rhs, [0.3; 0.5; -0.1], [1; 0.8], wrappedOpts);
assert(wrappedResult.has_orbit, 'Expected wrapped RHS + matcont_odefile flow to succeed.');
assert(wrappedResult.success, 'Expected wrapped RHS + matcont_odefile flow to succeed.');
assert(abs(wrappedResult.period - directResult.period) < 0.5, ...
    'Expected wrapped and direct MATCONT flows to agree on the period.');
assert(abs(wrappedResult.output_active_parameter_value - wrappedResult.input_active_parameter_value) <= wrappedOpts.matcont_parameter_tolerance, ...
    'Expected wrapped flow to return the input active parameter within tolerance.');

didThrow = false;
try
    orbitmatcont.extract_periodic_orbit(rhs, [0.3; 0.5; -0.1], [1; 0.8], rmfield(wrappedOpts, 'matcont_odefile'));
catch ME
    didThrow = strcmp(ME.identifier, 'orbitmatcont:MissingMatcontOdefile');
end
assert(didThrow, 'Expected a MissingMatcontOdefile error when matcont_odefile is omitted for a plain RHS.');

badNcolOpts = wrappedOpts;
badNcolOpts.matcont_ncol = 8;
didThrow = false;
try
    orbitmatcont.extract_periodic_orbit(rhs, [0.3; 0.5; -0.1], [1; 0.8], badNcolOpts);
catch ME
    didThrow = strcmp(ME.identifier, 'orbitmatcont:InvalidNcol');
end
assert(didThrow, 'Expected an InvalidNcol error when ncol falls outside the MATCONT-supported range.');

strictOpts = wrappedOpts;
strictOpts.matcont_parameter_tolerance = 1e-6;
strictResult = orbitmatcont.extract_periodic_orbit(rhs, [0.3; 0.5; -0.1], [1; 0.8], strictOpts);
assert(strictResult.success, 'Expected userfunction-based target continuation to recover the input parameter at 1e-6 tolerance.');
assert(abs(strictResult.output_active_parameter_value - strictResult.input_active_parameter_value) <= strictOpts.matcont_parameter_tolerance, ...
    'Expected strictResult to satisfy the tighter active-parameter tolerance.');

failureOpts = wrappedOpts;
failureOpts.matcont_parameter_tolerance = 1e-12;
failureOpts.matcont_return_max_points = 2;
failureOpts.matcont_return_scan_both_directions = false;
didThrow = false;
try
    orbitmatcont.extract_periodic_orbit(rhs, [0.3; 0.5; -0.1], [1; 0.8], failureOpts);
catch ME
    didThrow = strcmp(ME.identifier, 'orbitmatcont:ParameterReturnFailed');
end
assert(didThrow, 'Expected a ParameterReturnFailed error when the target continuation budget is too short to reach the userfunction zero.');

disp('All orbitmatcont tests passed.');
end
