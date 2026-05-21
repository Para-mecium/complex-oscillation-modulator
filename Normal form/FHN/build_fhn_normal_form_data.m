function data = build_fhn_normal_form_data(opts)
%BUILD_FHN_NORMAL_FORM_DATA Generate FHN AM/FM normal-form curves without plotting.

ensure_paths();

if nargin < 1
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);
core = build_fhn_core(cfg);
model = normalform.build_symbolic_two_state_model(build_fhn_spec(cfg, core));
core.auto_omega0 = normalform.frequency_from_linearization(model, 0, 0);
core.auto_chi0 = normalform.first_lyapunov_from_tensors(model, 0, 0);

curves = solve_fhn_closed_form_curves(core, cfg);

checks = struct();
checks.am_omega_max_abs = max(abs(closed_form_omega_sq(core, curves.am.f11, curves.am.f12) - core.omega0_sq_closed), [], 'omitnan');
checks.fm_chi_max_abs = max(abs(closed_form_chi(core, curves.fm.f11, curves.fm.f12) - core.chi0), [], 'omitnan');
checks.constraint_max_abs = max([ ...
    max(abs(curves.am.f22 + curves.am.f11), [], 'omitnan'), ...
    max(abs(curves.fm.f22 + curves.fm.f11), [], 'omitnan'), ...
    max(abs(curves.am.f21 - core.a21_over_a12 .* curves.am.f12), [], 'omitnan'), ...
    max(abs(curves.fm.f21 - core.a21_over_a12 .* curves.fm.f12), [], 'omitnan')], [], 'omitnan');
checks = add_fhn_regression_checks(checks, model, core, curves, cfg);

data = struct();
data.meta = struct( ...
    'model', 'FHN', ...
    'parameters', struct( ...
        'theta', cfg.theta, ...
        'gamma', cfg.gamma, ...
        'I', cfg.I, ...
        'epsilon', cfg.epsilon), ...
    'grids', struct( ...
        'am_f11_range', cfg.am_f11_range(:), ...
        'fm_f11_range', cfg.fm_f11_range(:)), ...
    'solver', struct( ...
        'fm_search_interval', cfg.fm_search_interval), ...
    'symbolic_cache_enabled', true);
data.core = core;
data.curves = struct('am', curves.am, 'fm', curves.fm);
data.checks = checks;

if cfg.assert_match
    validateattributes(checks.am_omega_max_abs, {'double'}, {'<=', cfg.am_tolerance}, ...
        mfilename, 'AM omega residual');
    validateattributes(checks.fm_chi_max_abs, {'double'}, {'<=', cfg.fm_tolerance}, ...
        mfilename, 'FM chi residual');
    validateattributes(checks.constraint_max_abs, {'double'}, {'<=', cfg.constraint_tolerance}, ...
        mfilename, 'constraint residual');
end
end

function opts = default_options()
opts = struct();
opts.theta = 0.2;
opts.gamma = 2.5;
opts.I = 0.1;
opts.epsilon = 0.1;
opts.am_f11_range = (-0.4:0.01:0.4);
opts.fm_f11_range = (-0.2:0.01:0.4);
opts.fm_search_interval = [-4, 2];
opts.assert_match = true;
opts.am_tolerance = 1e-10;
opts.fm_tolerance = 1e-6;
opts.constraint_tolerance = 1e-12;
end

function core = build_fhn_core(cfg)
core = struct();
core.theta = cfg.theta;
core.gamma = cfg.gamma;
core.I = cfg.I;
core.epsilon = cfg.epsilon;
core.V0 = solve_equilibrium(core.theta, core.gamma, core.I);
core.W0 = core.V0 / core.gamma;
core.alpha = -core.theta + 2 * (core.theta + 1) * core.V0 - 3 * core.V0.^2;
core.beta = 1 + core.theta - 3 * core.V0;
core.mu = 0.5 * (core.alpha - core.epsilon * core.gamma);
core.omega0_sq_closed = closed_form_omega_sq(core, 0, 0);
core.omega0 = sqrt(core.omega0_sq_closed);
core.a12 = -1;
core.a21 = core.epsilon;
core.a21_over_a12 = core.a21 / core.a12;
core.jacobian = [core.alpha, core.a12; core.a21, -core.epsilon * core.gamma];
core.chi0 = closed_form_chi(core, 0, 0);
core.auto_omega0 = [];
core.auto_chi0 = [];
end

function spec = build_fhn_spec(cfg, core)
spec = struct();
spec.cache_key = sprintf('FHN|%.16g|%.16g|%.16g|%.16g|%.16g|%.16g', ...
    cfg.theta, cfg.gamma, cfg.I, cfg.epsilon, core.V0, core.W0);
spec.rhs_builder = @(u1, u2, f11, f12) fhn_shifted_rhs(u1, u2, f11, f12, cfg, core);
end

function rhs = fhn_shifted_rhs(u1, u2, f11, f12, cfg, core)
V = core.V0 + u1;
W = core.W0 + u2;
f21 = core.a21_over_a12 * f12;
f22 = -f11;

rhs = [ ...
    V .* (V - cfg.theta) .* (1 - V) - W + cfg.I + f11 .* u1 + f12 .* u2; ...
    cfg.epsilon .* (V - cfg.gamma .* W) + f21 .* u1 + f22 .* u2];
end

function V0 = solve_equilibrium(theta, gamma, I)
equilibriumEquation = @(V) gamma .* V .* (V - theta) .* (1 - V) - V + gamma .* I;
V0 = fzero(equilibriumEquation, 0.3);
end

function curves = solve_fhn_closed_form_curves(core, cfg)
amF11 = cfg.am_f11_range(:);
fmF11 = cfg.fm_f11_range(:);

amF12 = 1 - sqrt((core.omega0_sq_closed + xi_value(core, amF11).^2) ./ core.epsilon);
fmResidualFn = @(f11, f12) closed_form_chi(core, f11, f12) - core.chi0;
fmF12 = normalform.solve_curve_continuation(fmF11, fmResidualFn, 0, ...
    struct('search_interval', cfg.fm_search_interval));

curves = struct();
curves.am = normalform.make_control_table(amF11, amF12, core.a21_over_a12);
curves.fm = normalform.make_control_table(fmF11, fmF12, core.a21_over_a12);
end

function checks = add_fhn_regression_checks(checks, model, core, curves, cfg)
autoAmOmega = normalform.frequency_from_linearization(model, curves.am.f11, curves.am.f12);
autoFmChi = normalform.first_lyapunov_from_tensors(model, curves.fm.f11, curves.fm.f12);

closedAmOmega = sqrt(closed_form_omega_sq(core, curves.am.f11, curves.am.f12));
closedFmChi = closed_form_chi(core, curves.fm.f11, curves.fm.f12);

checks.am_closed_form_max_abs = max(abs(autoAmOmega - closedAmOmega), [], 'omitnan');
checks.fm_closed_form_max_abs = max(abs(autoFmChi - closedFmChi), [], 'omitnan');
checks.auto_core_omega_abs = abs(core.auto_omega0 - core.omega0);
checks.auto_core_chi_abs = abs(core.auto_chi0 - core.chi0);

legacyAmF12 = legacy_am_curve(curves.am.f11);
legacyFmF12 = legacy_fm_curve(curves.fm.f11, cfg.fm_search_interval);

checks.am_legacy_curve_max_abs = max(abs(curves.am.f12 - legacyAmF12), [], 'omitnan');
checks.fm_legacy_curve_max_abs = max(abs(curves.fm.f12 - legacyFmF12), [], 'omitnan');
end

function omegaSq = closed_form_omega_sq(core, f11, f12)
omegaSq = core.epsilon .* (1 - f12).^2 - xi_value(core, f11).^2;
end

function chi = closed_form_chi(core, f11, f12)
omegaSq = closed_form_omega_sq(core, f11, f12);
denom1 = core.mu.^2 + omegaSq;
denom2 = core.mu.^2 + 9 * omegaSq;
chi = -3/8 + (core.beta.^2 ./ (4 * denom1)) .* ...
    (f11 + 2 * core.alpha - core.epsilon * core.gamma - ...
    4 * core.mu .* (core.epsilon .* (1 - f12).^2 ./ denom2));
end

function xi = xi_value(core, f11)
xi = f11 + 0.5 * (core.alpha + core.epsilon * core.gamma);
end

function f12 = legacy_am_curve(f11)
radicand = (32 * f11.^2) / 5 + ...
    (2269549760824611 .* f11) / 703687441776640 + ...
    10141204802429623294745511139173 / 15845632502852867518708790067200;
f12 = 1 - (5 / 4) .* sqrt(radicand);
end

function f12 = legacy_fm_curve(f11Grid, searchInterval)
residualFn = @(f11, f12) legacy_fm_residual(f11, f12);
f12 = normalform.solve_curve_continuation(f11Grid, residualFn, 0, ...
    struct('search_interval', searchInterval));
end

function value = legacy_fm_residual(f11, f12)
numerator = ...
    0.004499999999 .* f12.^7 - 0.03149999999 .* f12.^6 + ...
    (-0.1350000000 .* f11.^2 - 0.07037624749 .* f11 + 0.08532676680) .* f12.^5 + ...
    (0.6749999999 .* f11.^2 + 0.3518812377 .* f11 - 0.1116338340) .* f12.^4 + ...
    (1.350000000 .* f11.^4 + 1.40752495000000 .* f11.^3 - 0.608320764900000 .* f11 - ...
    0.800047835100000 .* f11.^2 + 0.0719752263400000) .* f12.^3 + ...
    (-4.049999999 .* f11.^4 + 0.417437344500000 .* f11 - 0.299856494500000 .* f11.^2 - ...
    4.22257485100000 .* f11.^3 - 0.0213903430800000) .* f12.^2 + ...
    (-4.499999999 .* f11.^6 - 0.532198328200001 .* f11.^4 - 7.03762474900000 .* f11.^5 - ...
    0.0977896198799999 .* f11 + 0.664759008600000 .* f11.^2 + 2.63262176000000 .* f11.^3 + ...
    0.00286137125000000) .* f12 + ...
    0.182428140399999 .* f11.^3 + 3.23219832900000 .* f11.^4 + 7.03762474900000 .* f11.^5 + ...
    4.499999999 .* f11.^6 - 0.000139187324600000 - 0.104854679100000 .* f11.^2 + ...
    0.00716805035899999 .* f11;

denominator = ...
    (-0.1095523154 + 1.511749833 .* f11 + 3 .* f11.^2 + 0.6 .* f12 - 0.3 .* f12.^2) .* ...
    (-0.03652084731 + 0.503916611 .* f11 + f11.^2 + 0.2 .* f12 - 0.1 .* f12.^2) .* ...
    (-1 + f12) .* ...
    (-0.1460680494 + 2.015666444 .* f11 + 4 .* f11.^2 + 0.8 .* f12 - 0.4 .* f12.^2);

value = numerator ./ denominator + 0.2381672552;
end
