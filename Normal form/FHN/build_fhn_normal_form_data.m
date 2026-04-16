function data = build_fhn_normal_form_data(opts)
%BUILD_FHN_NORMAL_FORM_DATA Generate FHN AM/FM normal-form curves without plotting.

ensure_paths();

if nargin < 1
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);
core = build_fhn_core(cfg);

amF11 = cfg.am_f11_range(:);
fmF11 = cfg.fm_f11_range(:);

amF12 = am_curve_from_normal_form(core, amF11);
fmResidualFn = @(f11, f12) chi_from_normal_form(core, f11, f12) - core.chi0;
fmF12 = normalform.solve_curve_continuation(fmF11, fmResidualFn, 0, ...
    struct('search_interval', cfg.fm_search_interval));

amTable = normalform.make_control_table(amF11, amF12, core.a21_over_a12);
fmTable = normalform.make_control_table(fmF11, fmF12, core.a21_over_a12);

checks = struct();
checks.am_omega_max_abs = max(abs(omega_sq(core, amTable.f11, amTable.f12) - core.omega0_sq), [], 'omitnan');
checks.fm_chi_max_abs = max(abs(fmResidualFn(fmTable.f11, fmTable.f12)), [], 'omitnan');
checks.constraint_max_abs = max([ ...
    max(abs(amTable.f22 + amTable.f11), [], 'omitnan'), ...
    max(abs(fmTable.f22 + fmTable.f11), [], 'omitnan'), ...
    max(abs(amTable.f21 - core.a21_over_a12 .* amTable.f12), [], 'omitnan'), ...
    max(abs(fmTable.f21 - core.a21_over_a12 .* fmTable.f12), [], 'omitnan')], [], 'omitnan');

data = struct();
data.meta = struct( ...
    'model', 'FHN', ...
    'parameters', struct( ...
        'theta', cfg.theta, ...
        'gamma', cfg.gamma, ...
        'I', cfg.I, ...
        'epsilon', cfg.epsilon), ...
    'grids', struct( ...
        'am_f11_range', amF11, ...
        'fm_f11_range', fmF11), ...
    'solver', struct( ...
        'fm_search_interval', cfg.fm_search_interval));
data.core = core;
data.curves = struct('am', amTable, 'fm', fmTable);
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
opts.am_f11_range = (-0.5:0.01:0.8);
opts.fm_f11_range = (-0.25:0.01:0.8);
opts.fm_search_interval = [-4, 2];
opts.assert_match = true;
opts.am_tolerance = 1e-10;
opts.fm_tolerance = 1e-6;
opts.constraint_tolerance = 1e-12;
end

function core = build_fhn_core(opts)
core = struct();
core.theta = opts.theta;
core.gamma = opts.gamma;
core.I = opts.I;
core.epsilon = opts.epsilon;
core.V0 = solve_equilibrium(core.theta, core.gamma, core.I);
core.W0 = core.V0 / core.gamma;
core.alpha = -core.theta + 2 * (core.theta + 1) * core.V0 - 3 * core.V0.^2;
core.beta = 1 + core.theta - 3 * core.V0;
core.mu = 0.5 * (core.alpha - core.epsilon * core.gamma);
core.omega0_sq = core.epsilon - ((core.alpha + core.epsilon * core.gamma) / 2).^2;
core.omega0 = sqrt(core.omega0_sq);
core.chi0 = chi_from_normal_form(core, 0, 0);
core.a12 = -1;
core.a21 = core.epsilon;
core.a21_over_a12 = core.a21 / core.a12;
end

function V0 = solve_equilibrium(theta, gamma, I)
equilibriumEquation = @(V) gamma .* V .* (V - theta) .* (1 - V) - V + gamma .* I;
V0 = fzero(equilibriumEquation, 0.3);
end

function f12 = am_curve_from_normal_form(core, f11)
xi = xi_value(core, f11);
f12 = 1 - sqrt((core.omega0_sq + xi.^2) ./ core.epsilon);
end

function omegaSq = omega_sq(core, f11, f12)
omegaSq = core.epsilon .* (1 - f12).^2 - xi_value(core, f11).^2;
end

function chi = chi_from_normal_form(core, f11, f12)
omegaSq = omega_sq(core, f11, f12);
denom1 = core.mu.^2 + omegaSq;
denom2 = core.mu.^2 + 9 * omegaSq;
chi = -3/8 + (core.beta.^2 ./ (4 * denom1)) .* ...
    (f11 + 2 * core.alpha - core.epsilon * core.gamma - ...
    4 * core.mu .* (core.epsilon .* (1 - f12).^2 ./ denom2));
end

function xi = xi_value(core, f11)
xi = f11 + 0.5 * (core.alpha + core.epsilon * core.gamma);
end
