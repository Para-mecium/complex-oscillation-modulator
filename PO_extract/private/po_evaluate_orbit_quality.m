function quality = po_evaluate_orbit_quality(result, rhsFun)
quality = struct( ...
    'closure_defect', NaN, ...
    'ode_residual_inf', NaN, ...
    'ode_residual_l2', NaN, ...
    'extrema_refinement_delta', NaN, ...
    'time_sample_count', 0, ...
    'poincare_residual_last', NaN, ...
    'parameter_error', NaN);

if isfield(result, 'diagnostics') && isstruct(result.diagnostics) && ...
        isfield(result.diagnostics, 'poincareResidual') && ~isempty(result.diagnostics.poincareResidual)
    quality.poincare_residual_last = result.diagnostics.poincareResidual(end);
end
if isfield(result, 'parameter_error') && isscalar(result.parameter_error) && isnumeric(result.parameter_error)
    quality.parameter_error = result.parameter_error;
end
if ~(isfield(result, 'has_orbit') && result.has_orbit)
    return;
end

t = result.orbit_t(:);
y = double(result.orbit_y);
quality.time_sample_count = size(y, 1);
if isempty(t) || isempty(y) || size(y, 1) ~= numel(t)
    return;
end

scale = max(1, max(abs(y(:))));
quality.closure_defect = norm(y(end, :) - y(1, :), inf) / scale;
quality.extrema_refinement_delta = extrema_refinement_delta(t, y);

if isa(rhsFun, 'function_handle') && numel(t) >= 3
    try
        rhs = zeros(size(y));
        for i = 1:numel(t)
            rhs(i, :) = reshape(rhsFun(t(i), y(i, :).'), 1, []);
        end
        dydt = estimate_derivative(t, y);
        residual = dydt - rhs;
        rowResidual = sqrt(sum(residual.^2, 2));
        quality.ode_residual_inf = max(rowResidual);
        quality.ode_residual_l2 = sqrt(mean(rowResidual.^2));
    catch
        quality.ode_residual_inf = NaN;
        quality.ode_residual_l2 = NaN;
    end
end

end

function dydt = estimate_derivative(t, y)
n = numel(t);
dydt = zeros(size(y));
if n < 2
    return;
end

dydt(1, :) = (y(2, :) - y(1, :)) / (t(2) - t(1));
dydt(end, :) = (y(end, :) - y(end - 1, :)) / (t(end) - t(end - 1));
for i = 2:(n - 1)
    dydt(i, :) = (y(i + 1, :) - y(i - 1, :)) / (t(i + 1) - t(i - 1));
end
end

function delta = extrema_refinement_delta(t, y)
delta = NaN;
if numel(t) < 3 || any(diff(t) <= 0)
    return;
end

rawMax = max(y, [], 1);
rawMin = min(y, [], 1);
denseCount = max(400, 5 * numel(t));
tq = linspace(t(1), t(end), denseCount).';
yq = interp1(t, y, tq, 'pchip');
denseMax = max(yq, [], 1);
denseMin = min(yq, [], 1);
delta = max(abs([denseMax - rawMax, denseMin - rawMin]), [], 'all');
end
