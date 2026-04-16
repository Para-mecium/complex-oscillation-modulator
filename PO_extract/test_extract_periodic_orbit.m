function test_extract_periodic_orbit()
%TEST_EXTRACT_PERIODIC_ORBIT Regression checks for PO_extract.

po_dir = fileparts(mfilename('fullpath'));
addpath(po_dir);

opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 20, ...
    'max_windows', 5);
stableFields = {'success', 'has_orbit', 'status', 'code', 'message', ...
    'period', 'orbit_t', 'orbit_y', 'amplitude', 'max_variable', 'min_variable'};

vdp_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, opts);
assert(all(isfield(vdp_result, stableFields)), 'Stable contract fields must always exist.');
assert(vdp_result.code == 2, 'Van der Pol oscillator should converge to a periodic orbit.');
assert(vdp_result.status == "converged_periodic_orbit", ...
    'Converged periodic orbit should report the stable status string.');
assert(vdp_result.success, 'Converged periodic orbit should set success=true.');
assert(vdp_result.has_orbit, 'Converged periodic orbit should expose an extracted orbit.');
assert(vdp_result.period > 5 && vdp_result.period < 8, 'Unexpected Van der Pol period.');
assert(size(vdp_result.orbit_y, 2) == 2, 'Orbit state dimension mismatch.');
assert(abs(vdp_result.orbit_t(end) - vdp_result.period) < 1e-10, 'Relative orbit time should end at the period.');

no_cycle_result = extract_periodic_orbit(@stable_linear_rhs, [1; -2], [], opts);
assert(any(no_cycle_result.code == [-1, 0]), 'Stable linear system should not produce a periodic orbit.');
assert(~no_cycle_result.success, 'Stable linear system should not set success=true.');
assert(~no_cycle_result.has_orbit, 'Stable linear system should not return an orbit.');
assert(isempty(no_cycle_result.period), 'No-orbit results should not expose a usable period.');
assert(isempty(no_cycle_result.orbit_t), 'No-orbit results should not expose a usable time grid.');
assert(isempty(no_cycle_result.amplitude), 'No-orbit results should not expose usable amplitudes.');
assert(isempty(no_cycle_result.max_variable) && isempty(no_cycle_result.min_variable), ...
    'No-orbit results should not expose usable extrema.');

short_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 3, ...
    'max_windows', 1);
too_short_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, short_opts);
assert(too_short_result.status == "no_periodic_orbit_detected_on_tspan", ...
    'Short horizon should map to the stable no-orbit status.');
assert(too_short_result.code == 0, 'Short horizon should not provide convincing periodic-orbit evidence.');
assert(~too_short_result.has_orbit, 'Short horizon should not return an orbit.');

drift_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 15, ...
    'max_windows', 6);
drift_result = extract_periodic_orbit(@oscillator_with_drift_rhs, [1; 0; 0], [], drift_opts);
assert(drift_result.status == "no_periodic_orbit_detected_on_tspan", ...
    'Drifting oscillator should report the stable no-orbit status.');
assert(drift_result.code == 0, 'Drifting oscillator should not be classified as a periodic orbit.');
assert(~drift_result.success, 'Drifting oscillator should not set success=true.');
assert(~drift_result.has_orbit, 'Drifting oscillator should not expose an orbit.');

candidate_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 14, ...
    'max_windows', 5);
candidate_result = extract_periodic_orbit(@slow_limit_cycle_rhs, [2.5; 0], [0.02; 1], candidate_opts);
assert(candidate_result.code == 1, 'Slowly converging limit cycle should produce a candidate orbit.');
assert(candidate_result.status == "candidate_periodic_orbit_not_converged", ...
    'Candidate periodic orbit should report the stable status string.');
assert(~candidate_result.success, 'Candidate orbit should not set success=true.');
assert(candidate_result.has_orbit, 'Candidate orbit should still return the last cycle.');
assert(~isempty(candidate_result.orbit_y), 'Candidate orbit data should be populated.');

decay_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 10, ...
    'max_windows', 3);
decay_result = extract_periodic_orbit(@stable_linear_rhs, [1e-8; -1e-8], [], decay_opts);
assert(decay_result.status == "decaying_to_equilibrium_or_nonoscillatory", ...
    'Small decaying trajectories should report the nonoscillatory stable status.');
assert(decay_result.code == -1, 'Small decaying trajectories should map to code -1.');
assert(~decay_result.has_orbit, 'Nonoscillatory trajectories should not expose an orbit.');

negative_result = extract_periodic_orbit(@shifted_oscillator_rhs, [-1; 0], 0, drift_opts);
assert(negative_result.code == 2, 'Shifted oscillator should still converge to a periodic orbit.');
assert(negative_result.max_variable(1) < 0, 'Shifted oscillator should keep the first-state maxima negative.');

row_vector_result = extract_periodic_orbit(@vdp_rhs, [2, 0], 1, opts);
assert(row_vector_result.code == 2, 'Row-vector y0 should be normalized automatically.');

bad_solver_opts = opts;
bad_solver_opts.solver_name = 'ode_not_available';
bad_solver_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, bad_solver_opts);
assert(bad_solver_result.status == "invalid_input", 'Invalid solver names must fail validation.');
assert(isnan(bad_solver_result.code), 'Invalid input should leave the numeric code undefined.');
assert(~bad_solver_result.success && ~bad_solver_result.has_orbit, ...
    'Invalid input should not report success or orbit data.');

solver_option_opts = opts;
solver_option_opts.solver_options = odeset('Events', @premature_stop_event);
solver_option_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, solver_option_opts);
assert(solver_option_result.code == 2, ...
    'Internal event detection should override external solver_options.Events.');

bad_solver_option_opts = opts;
bad_solver_option_opts.solver_options = 7;
bad_solver_option_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, bad_solver_option_opts);
assert(bad_solver_option_result.status == "invalid_input", ...
    'Non-struct solver_options must fail validation.');

invalid_option_opts = opts;
invalid_option_opts.event = "bad_event_option";
invalid_option_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, invalid_option_opts);
assert(invalid_option_result.status == "invalid_options", ...
    'Malformed legacy event options must report invalid_options.');
assert(isnan(invalid_option_result.code), 'invalid_options should leave the numeric code undefined.');
assert(~invalid_option_result.success && ~invalid_option_result.has_orbit, ...
    'invalid_options should not report success or orbit data.');

solver_failure_opts = opts;
solver_failure_opts.solver = @failing_solver;
solver_failure_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, solver_failure_opts);
assert(solver_failure_result.status == "solver_failed", ...
    'Solver execution failures must report solver_failed.');
assert(isnan(solver_failure_result.code), 'solver_failed should leave the numeric code undefined.');
assert(~solver_failure_result.success && ~solver_failure_result.has_orbit, ...
    'solver_failed should not report success or orbit data.');

detection_failure_opts = opts;
detection_failure_opts.observableFcn = @bad_observable;
detection_failure_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, detection_failure_opts);
assert(detection_failure_result.status == "detection_failed", ...
    'Observable failures outside the solver should report detection_failed.');
assert(isnan(detection_failure_result.code), 'detection_failed should leave the numeric code undefined.');
assert(~detection_failure_result.success && ~detection_failure_result.has_orbit, ...
    'detection_failed should not report success or orbit data.');

legacy_event_opts = opts;
legacy_event_opts.event = 2;
legacy_event_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, legacy_event_opts);
assert(legacy_event_result.code == 2 && legacy_event_result.has_orbit, ...
    'Legacy numeric event selection should remain supported.');

legacy_min_events_opts = opts;
legacy_min_events_opts.min_events = 9;
legacy_min_events_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, legacy_min_events_opts);
assert(legacy_min_events_result.diagnostics.thresholds.minCrossings == 9, ...
    'Legacy min_events should still map to minCrossings.');

legacy_match_tol_opts = opts;
legacy_match_tol_opts.match_tol = 1e-3;
legacy_match_tol_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, legacy_match_tol_opts);
assert(abs(legacy_match_tol_result.diagnostics.thresholds.strict.poincareTol - 1e-3) < 1e-12, ...
    'Legacy match_tol should still map to poincareTol.');

legacy_period_tol_opts = opts;
legacy_period_tol_opts.period_tol = 2e-3;
legacy_period_tol_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, legacy_period_tol_opts);
assert(abs(legacy_period_tol_result.diagnostics.thresholds.strict.periodTol - 2e-3) < 1e-12, ...
    'Legacy period_tol should still map to periodTol.');

matcont_opts = drift_opts;
matcont_opts.backend = 'matcont';
matcont_opts.matcont_root = fullfile(fileparts(po_dir), 'MatCont7p6');
matcont_result = extract_periodic_orbit(@shifted_oscillator_rhs, [-1; 0], 0, matcont_opts);
assert(matcont_result.code == 2, 'Direct detection should still work when MATCONT is requested.');
assert(matcont_result.backend_used == "direct", 'Backend should remain direct when MATCONT is disabled.');
assert(isfield(matcont_result.diagnostics, 'matcont') && ...
    matcont_result.diagnostics.matcont.status == "disabled", ...
    'MATCONT diagnostics should report disabled status.');
assert(matcont_result.diagnostics.matcont.requested, ...
    'MATCONT diagnostics should record that refinement was requested.');

api_opts = struct('solver_name', 'ode45', 'tspan', [0, 100]);
api_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, api_opts);
assert(api_result.code == 2, 'extract_periodic_orbit should return a converged periodic orbit when available.');
assert(api_result.status == "converged_periodic_orbit", ...
    'Converged orbit should report the new structured status.');
assert(api_result.success, 'Converged orbit should set success=true.');
assert(api_result.has_orbit, 'Converged orbit should expose orbit data.');
assert(~isempty(api_result.orbit_t) && ~isempty(api_result.orbit_y), ...
    'Converged orbit should populate the extracted orbit arrays.');
assert(api_result.period > 5 && api_result.period < 8, ...
    'extract_periodic_orbit should report the expected Van der Pol period.');
assert(all(api_result.amplitude > 0), 'Structured amplitude should be positive on periodic data.');
assert(all(api_result.max_variable > api_result.min_variable), ...
    'Structured extrema should be ordered.');

candidate_api_opts = struct('solver_name', 'ode45', 'tspan', [0, 70]);
candidate_api_result = extract_periodic_orbit(@slow_limit_cycle_rhs, [2.5; 0], [0.02; 1], candidate_api_opts);
assert(candidate_api_result.code == 1, ...
    'Slowly converging limit cycle should report a candidate periodic orbit.');
assert(candidate_api_result.status == "candidate_periodic_orbit_not_converged", ...
    'Candidate orbit should report the new structured status.');
assert(~candidate_api_result.success, 'Candidate orbit should keep success=false.');
assert(candidate_api_result.has_orbit, 'Candidate orbit should still expose orbit data.');
assert(~isempty(candidate_api_result.orbit_t) && ~isempty(candidate_api_result.orbit_y), ...
    'Candidate orbit should populate the extracted orbit arrays.');
assert(~isempty(candidate_api_result.amplitude) && ~isempty(candidate_api_result.period), ...
    'Candidate orbit should return structured observables together with the orbit.');

disp('All PO_extract tests passed.');
end

function dydt = vdp_rhs(~, y, mu)
y = y(:);
dydt = [y(2); mu * (1 - y(1)^2) * y(2) - y(1)];
end

function dydt = stable_linear_rhs(~, y, ~)
y = y(:);
dydt = [-y(1); -2 * y(2)];
end

function dydt = oscillator_with_drift_rhs(~, y, ~)
y = y(:);
dydt = [y(2); -y(1); 0.5];
end

function dydt = slow_limit_cycle_rhs(~, y, parameter)
y = y(:);
alpha = parameter(1);
omega = parameter(2);
r = hypot(y(1), y(2));
growth = alpha * (1 - r);
dydt = [growth * y(1) - omega * y(2); growth * y(2) + omega * y(1)];
end

function dydt = shifted_oscillator_rhs(~, y, omega)
y = y(:);
z1 = y(1) + 2;
dydt = [y(2); -(omega + 1)^2 * z1];
end

function [value, isterminal, direction] = premature_stop_event(t, ~)
value = t - 0.1;
isterminal = 1;
direction = 1;
end

function varargout = failing_solver(varargin)
error('test:FailingSolver', 'Intentional solver failure for contract regression.');
end

function value = bad_observable(~, ~, ~, ~)
value = [1; 2];
end
