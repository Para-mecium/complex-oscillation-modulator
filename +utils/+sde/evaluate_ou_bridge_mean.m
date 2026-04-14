function xi = evaluate_ou_bridge_mean(t, tGrid, xiGrid, tau)
tGrid = reshape(double(tGrid), [], 1);

if nargin < 4 || isempty(tau)
    tau = 1;
end
tau = double(tau);
if ~(isscalar(tau) && tau > 0)
    error('utils:sde:evaluate_ou_bridge_mean:InvalidTau', ...
        'tau must be a positive scalar.');
end
if size(xiGrid, 1) ~= numel(tGrid)
    error('utils:sde:evaluate_ou_bridge_mean:SizeMismatch', ...
        'xiGrid must have numel(tGrid) rows.');
end
if isempty(tGrid)
    error('utils:sde:evaluate_ou_bridge_mean:EmptyGrid', ...
        'tGrid must contain at least one point.');
end

if t <= tGrid(1)
    xi = xiGrid(1, :).';
    return
end
if t >= tGrid(end)
    xi = xiGrid(end, :).';
    return
end

idx = find(tGrid <= t, 1, 'last');
idx = min(max(idx, 1), numel(tGrid) - 1);

t0 = tGrid(idx);
t1 = tGrid(idx + 1);
h = t1 - t0;
s = t - t0;
lambda = 1 / tau;

if lambda * h < 1e-6
    w0 = 1 - s / h;
    w1 = s / h;
else
    denom = sinh(lambda * h);
    w0 = sinh(lambda * (h - s)) / denom;
    w1 = sinh(lambda * s) / denom;
end

xi = (w0 * xiGrid(idx, :) + w1 * xiGrid(idx + 1, :)).';
end
