function xiGrid = generate_ou_path(tGrid, sigma, tau)
tGrid = reshape(double(tGrid), [], 1);
sigma = reshape(double(sigma), 1, []);

if isempty(tGrid)
    error('utils:sde:generate_ou_path:EmptyGrid', ...
        'tGrid must contain at least one point.');
end
if nargin < 3 || isempty(tau)
    tau = 1;
end
tau = double(tau);
if ~(isscalar(tau) && tau > 0)
    error('utils:sde:generate_ou_path:InvalidTau', ...
        'tau must be a positive scalar.');
end
if any(diff(tGrid) < 0)
    error('utils:sde:generate_ou_path:NonMonotoneGrid', ...
        'tGrid must be nondecreasing.');
end

nSteps = numel(tGrid) - 1;
nNoise = numel(sigma);
xiGrid = zeros(numel(tGrid), nNoise);
if nSteps == 0
    return
end

dtVec = diff(tGrid);
for k = 1:nSteps
    alpha = exp(-dtVec(k) / tau);
    noiseStd = sigma .* sqrt((1 - alpha^2) / (2 * tau));
    xiGrid(k + 1, :) = alpha * xiGrid(k, :) + noiseStd .* randn(1, nNoise);
end
end
