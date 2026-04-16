function curves = solve_am_fm_curves(model, core, opts)
%SOLVE_AM_FM_CURVES Solve AM/FM control curves for a two-state model.

if nargin < 3
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);

amF11 = cfg.am_f11_range(:);
fmF11 = cfg.fm_f11_range(:);

amResidualFn = @(f11, f12) normalform.frequency_from_linearization(model, f11, f12) - core.omega0;
fmResidualFn = @(f11, f12) normalform.first_lyapunov_from_tensors(model, f11, f12) - core.chi0;

solverOpts = struct('search_interval', cfg.search_interval);
amF12 = normalform.solve_curve_continuation(amF11, amResidualFn, cfg.am_seed, solverOpts);
fmF12 = normalform.solve_curve_continuation(fmF11, fmResidualFn, cfg.fm_seed, solverOpts);

amF12 = polish_curve_points(amF11, amF12, amResidualFn, cfg.search_interval);
fmF12 = polish_curve_points(fmF11, fmF12, fmResidualFn, cfg.search_interval);

curves = struct();
curves.am = normalform.make_control_table(amF11, amF12, core.a21_over_a12);
curves.fm = normalform.make_control_table(fmF11, fmF12, core.a21_over_a12);
end

function cfg = default_options()
cfg = struct();
cfg.am_f11_range = [];
cfg.fm_f11_range = [];
cfg.search_interval = [-4, 4];
cfg.am_seed = 0;
cfg.fm_seed = 0;
end

function curve = polish_curve_points(f11Grid, initialCurve, residualFn, searchInterval)
curve = initialCurve;
for idx = 1:numel(f11Grid)
    curve(idx) = polish_root(f11Grid(idx), curve(idx), residualFn, searchInterval);
end
end

function root = polish_root(f11, guess, residualFn, searchInterval)
if ~isfinite(guess)
    root = guess;
    return
end

if exist('fsolve', 'file') == 2
    solverOptions = optimoptions('fsolve', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-14, ...
        'StepTolerance', 1e-14, ...
        'OptimalityTolerance', 1e-14);
    root = fsolve(@(f12) residualFn(f11, f12), guess, solverOptions);
    return
end

root = fzero(@(f12) residualFn(f11, f12), ...
    bracket_near_guess(f11, guess, residualFn, searchInterval));
end

function bracket = bracket_near_guess(f11, guess, residualFn, searchInterval)
span = 1e-3;
for iter = 1:10
    left = max(searchInterval(1), guess - span);
    right = min(searchInterval(2), guess + span);
    leftValue = residualFn(f11, left);
    rightValue = residualFn(f11, right);
    if isfinite(leftValue) && isfinite(rightValue) && leftValue * rightValue <= 0
        bracket = [left, right];
        return
    end
    span = span * 2;
end

error('normalform:PolishBracketFailure', ...
    'Failed to bracket a root while polishing the continuation result.');
end
