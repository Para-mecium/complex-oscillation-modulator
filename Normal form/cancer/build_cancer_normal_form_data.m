function data = build_cancer_normal_form_data(opts)
%BUILD_CANCER_NORMAL_FORM_DATA Generate cancer-example AM/FM normal-form curves.

ensure_paths();

if nargin < 1
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);
core = build_cancer_core(cfg);
model = normalform.build_symbolic_two_state_model(build_cancer_spec(cfg, core));
core.chi0 = normalform.first_lyapunov_from_tensors(model, 0, 0);

curves = normalform.solve_am_fm_curves(model, core, struct( ...
    'am_f11_range', cfg.am_f11_range, ...
    'fm_f11_range', cfg.fm_f11_range, ...
    'search_interval', cfg.fm_search_interval, ...
    'am_seed', 0, ...
    'fm_seed', 0));

checks = normalform.compute_two_state_checks(model, core, curves.am, curves.fm);
checks.equilibrium_residual_max_abs = max(abs(uncontrolled_rhs(cfg, cfg.phi0, cfg.mu0)), [], 'omitnan');

data = struct();
data.meta = struct( ...
    'model', 'cancer', ...
    'parameters', struct( ...
        'alpha', cfg.alpha, ...
        'kappa', cfg.kappa, ...
        'gamma1', cfg.gamma1, ...
        'gamma2', cfg.gamma2, ...
        'epsilon', cfg.epsilon, ...
        'phi0', cfg.phi0, ...
        'mu0', cfg.mu0, ...
        'K', cfg.K), ...
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
    validateattributes(checks.equilibrium_residual_max_abs, {'double'}, ...
        {'<=', cfg.equilibrium_tolerance}, mfilename, 'equilibrium residual');
end
end

function opts = default_options()
opts = struct();
opts.alpha = 0.35;
opts.kappa = 5;
opts.gamma1 = 1;
opts.gamma2 = 2.5;
opts.epsilon = 0.08;
opts.phi0 = 7 / 4;
opts.mu0 = 11 / 4;
opts.K = 5;
opts.am_f11_range = (-2.0:0.01:3.0);
opts.fm_f11_range = (-0.5:0.01:2.0);
opts.fm_search_interval = [-8, 8];
opts.assert_match = true;
opts.am_tolerance = 1e-10;
opts.fm_tolerance = 1e-8;
opts.constraint_tolerance = 1e-12;
opts.equilibrium_tolerance = 1e-10;
end

function core = build_cancer_core(cfg)
[phiEq, muEq] = solve_cancer_equilibrium(cfg);

if max(abs([phiEq - cfg.phi0; muEq - cfg.mu0])) > 100 * cfg.equilibrium_tolerance
    error('build_cancer_normal_form_data:EquilibriumMismatch', ...
        ['Configured (phi0, mu0) does not match the equilibrium implied by ', ...
         'the model parameters.']);
end

jacobian0 = uncontrolled_jacobian(cfg, phiEq, muEq);
a12 = jacobian0(1, 2);
a21 = jacobian0(2, 1);

if abs(a12) <= eps
    error('build_cancer_normal_form_data:InvalidJacobian', ...
        'a12 is zero at the equilibrium, so the structural ratio is undefined.');
end

core = struct();
core.phi0 = phiEq;
core.mu0 = muEq;
core.K = cfg.K;
core.jacobian = jacobian0;
core.a12 = a12;
core.a21 = a21;
core.a21_over_a12 = a21 / a12;
core.omega0 = sqrt(det(jacobian0) - 0.25 * trace(jacobian0)^2);
core.chi0 = [];
end

function spec = build_cancer_spec(cfg, core)
spec = struct();
spec.cache_key = sprintf('cancer|%.16g|%.16g|%.16g|%.16g|%.16g|%.16g|%.16g|%.16g|%.16g', ...
    cfg.alpha, cfg.kappa, cfg.gamma1, cfg.gamma2, cfg.epsilon, ...
    core.phi0, core.mu0, cfg.K, core.a21_over_a12);
spec.rhs_builder = @(u1, u2, f11, f12) cancer_shifted_rhs(u1, u2, f11, f12, cfg, core);
end

function rhs = cancer_shifted_rhs(u1, u2, f11, f12, cfg, core)
phi = core.phi0 + u1;
mu = core.mu0 + u2;

mmPhi = u1 / (u1 + cfg.K);
mmMu = u2 / (u2 + cfg.K);
f21 = core.a21_over_a12 * f12;
f22 = -f11;

rhs = [ ...
    (cfg.alpha + cfg.kappa * phi.^2 ./ (cfg.gamma1 + phi.^2 + cfg.gamma2 * mu) - phi) / cfg.epsilon + ...
    f11 * mmPhi + f12 * mmMu; ...
    1 + phi - mu + f21 * mmPhi + f22 * mmMu];
end

function [phiEq, muEq] = solve_cancer_equilibrium(cfg)
steadyStateEquation = @(phi) cfg.alpha + ...
    cfg.kappa .* phi.^2 ./ (cfg.gamma1 + phi.^2 + cfg.gamma2 .* (1 + phi)) - phi;
phiEq = fzero(steadyStateEquation, cfg.phi0);
muEq = 1 + phiEq;
end

function rhs = uncontrolled_rhs(cfg, phi, mu)
rhs = [ ...
    (cfg.alpha + cfg.kappa .* phi.^2 ./ (cfg.gamma1 + phi.^2 + cfg.gamma2 .* mu) - phi) ./ cfg.epsilon; ...
    1 + phi - mu];
end

function jacobian0 = uncontrolled_jacobian(cfg, phi0, mu0)
denom = cfg.gamma1 + phi0.^2 + cfg.gamma2 .* mu0;
jacobian0 = [ ...
    (-1 + 2 * cfg.kappa * phi0 .* (cfg.gamma1 + cfg.gamma2 .* mu0) ./ denom.^2) ./ cfg.epsilon, ...
    -(cfg.kappa * cfg.gamma2 .* phi0.^2 ./ denom.^2) ./ cfg.epsilon; ...
    1, -1];
end
