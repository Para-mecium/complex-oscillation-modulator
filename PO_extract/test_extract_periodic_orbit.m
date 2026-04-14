function test_extract_periodic_orbit()
%TEST_EXTRACT_PERIODIC_ORBIT Regression checks for PO_extract.

po_dir = fileparts(mfilename('fullpath'));
addpath(po_dir);

opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 20, ...
    'max_windows', 5);

vdp_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, opts);
assert(vdp_result.code == 2, 'Van der Pol oscillator should converge to a periodic orbit.');
assert(vdp_result.success, 'Converged periodic orbit should set success=true.');
assert(vdp_result.has_orbit, 'Converged periodic orbit should expose an extracted orbit.');
assert(vdp_result.period > 5 && vdp_result.period < 8, 'Unexpected Van der Pol period.');
assert(size(vdp_result.orbit_y, 2) == 2, 'Orbit state dimension mismatch.');
assert(abs(vdp_result.orbit_t(end) - vdp_result.period) < 1e-10, 'Relative orbit time should end at the period.');

no_cycle_result = extract_periodic_orbit(@stable_linear_rhs, [1; -2], [], opts);
assert(any(no_cycle_result.code == [-1, 0]), 'Stable linear system should not produce a periodic orbit.');
assert(~no_cycle_result.success, 'Stable linear system should not set success=true.');
assert(~no_cycle_result.has_orbit, 'Stable linear system should not return an orbit.');

short_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 3, ...
    'max_windows', 1);
too_short_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, short_opts);
assert(too_short_result.code == 0, 'Short horizon should not provide convincing periodic-orbit evidence.');
assert(~too_short_result.has_orbit, 'Short horizon should not return an orbit.');

drift_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 15, ...
    'max_windows', 6);
drift_result = extract_periodic_orbit(@oscillator_with_drift_rhs, [1; 0; 0], [], drift_opts);
assert(drift_result.code == 0, 'Drifting oscillator should not be classified as a periodic orbit.');
assert(~drift_result.success, 'Drifting oscillator should not set success=true.');

candidate_opts = struct( ...
    'solver_name', 'ode45', ...
    'single_timespan', 14, ...
    'max_windows', 5);
candidate_result = extract_periodic_orbit(@slow_limit_cycle_rhs, [2.5; 0], [0.02; 1], candidate_opts);
assert(candidate_result.code == 1, 'Slowly converging limit cycle should produce a candidate orbit.');
assert(~candidate_result.success, 'Candidate orbit should not set success=true.');
assert(candidate_result.has_orbit, 'Candidate orbit should still return the last cycle.');
assert(~isempty(candidate_result.orbit_y), 'Candidate orbit data should be populated.');

negative_result = extract_periodic_orbit(@shifted_oscillator_rhs, [-1; 0], 0, drift_opts);
assert(negative_result.code == 2, 'Shifted oscillator should still converge to a periodic orbit.');
assert(negative_result.max_variable(1) < 0, 'Shifted oscillator should keep the first-state maxima negative.');

row_vector_result = extract_periodic_orbit(@vdp_rhs, [2, 0], 1, opts);
assert(row_vector_result.code == 2, 'Row-vector y0 should be normalized automatically.');

bad_solver_opts = opts;
bad_solver_opts.solver_name = 'ode_not_available';
bad_solver_result = extract_periodic_orbit(@vdp_rhs, [2; 0], 1, bad_solver_opts);
assert(bad_solver_result.status == "invalid_input", 'Invalid solver names must fail validation.');

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

[TS, amplitude_variable, period, max_variable, min_variable] = POinfo(@vdp_rhs, [2; 0], 1, 'ode45', 100);
assert(~isempty(TS{1, 1}) && ~isempty(TS{1, 2}), 'POinfo should return a converged periodic orbit when available.');
assert(period > 5 && period < 8, 'POinfo should report the expected Van der Pol period.');
assert(all(amplitude_variable > 0), 'POinfo amplitude should be positive on periodic data.');
assert(all(max_variable > min_variable), 'POinfo extrema should be ordered.');

lastwarn('', '');
[TS_candidate, amp_candidate, period_candidate] = POinfo(@slow_limit_cycle_rhs, [2.5; 0], [0.02; 1], 'ode45', 70);
[warnMsg, warnId] = lastwarn();
assert(~isempty(TS_candidate{1, 1}) && ~isempty(TS_candidate{1, 2}), ...
    'POinfo should return candidate orbit data when strict convergence fails.');
assert(~isempty(amp_candidate) && ~isempty(period_candidate), ...
    'POinfo should return candidate observables together with the orbit.');
assert(strcmp(warnId, 'ErgodicMethod_ODE:PeriodicOrbitNotFound'), ...
    'POinfo candidate path should still emit the legacy warning identifier.');
assert(~isempty(warnMsg), 'POinfo candidate path should emit a warning message.');

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
