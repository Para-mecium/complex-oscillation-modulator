function result = Circadian_SDE(params, y0, opts)
ensure_shared_sde_path();

if nargin < 3 || isempty(opts)
    opts = struct();
end
if numel(opts) > 1
    noiseClass = {opts.noiseClass};
    opts = opts(1);
    opts.noiseClass = noiseClass;
end

if ~isfield(opts, 'T')
    opts.T = 60;
end
if ~isfield(opts, 'dt')
    opts.dt = 0.01;
end
if ~isfield(opts, 'sigma')
    opts.sigma = 0.4;
end
if ~isfield(opts, 'seed')
    opts.seed = 1;
end
if ~isfield(opts, 'noiseClass')
    opts.noiseClass = repmat({'o'}, 1, 3);
end
if ~isfield(opts, 'solver')
    opts.solver = 'ou_ode45';
end
if ~isfield(opts, 'odeRelTol')
    opts.odeRelTol = 1e-8;
end
if ~isfield(opts, 'odeAbsTol')
    opts.odeAbsTol = 1e-10;
end

params = reshape(double(params), 1, []);
y0 = reshape(double(y0), [], 1);
nNoise = 3;
if numel(params) ~= 2
    error('Circadian_SDE:InvalidParams', ...
        'Expected params to be [Kd, AT].');
end
if numel(y0) ~= 3
    error('Circadian_SDE:InvalidInitialState', ...
        'Expected y0 to contain [M; Pc; Pn].');
end

N = round(opts.T / opts.dt);
sigmaVec = normalize_sigma(opts.sigma, nNoise);
rng(double(opts.seed), 'twister');

if all(sigmaVec == 0)
    [t, X] = simulate_deterministic(params, y0, opts.T, N, opts);
    solverUsed = 'ode45_deterministic';
    ouGridSize = [0, 0];
else
    switch char(opts.solver)
        case 'ou_ode45'
            [t, X, ouGridSize] = simulate_stochastic_with_ode45(params, y0, opts.T, N, sigmaVec, opts);
            solverUsed = 'ou_ode45';
        case 'milstein_ou'
            [t, X, ouGridSize] = simulate_stochastic_with_milstein(params, y0, opts.T, N, sigmaVec, opts);
            solverUsed = 'milstein_ou';
        otherwise
            error('Circadian_SDE:UnknownSolver', ...
                'Unknown solver "%s".', char(opts.solver));
    end
end

t = t(:);
if size(X, 1) ~= numel(t) || size(X, 2) ~= 3
    error('Circadian_SDE:InvalidOutputShape', ...
        'Expected X to have size [numel(t), 3].');
end

result = struct();
result.t = t;
result.X = X;
result.params = params;
result.y0 = y0;
result.opts = opts;
result.meta = struct( ...
    'seed', double(opts.seed), ...
    'solver', solverUsed, ...
    'dt', double(opts.dt), ...
    'numSteps', N, ...
    'ouGridSize', reshape(double(ouGridSize), 1, []));
end

function [t, X] = simulate_deterministic(params, y0, totalTime, nSteps, opts)
t = linspace(0, totalTime, nSteps + 1).';
odeOpts = odeset('RelTol', opts.odeRelTol, 'AbsTol', opts.odeAbsTol);
[~, Y] = ode45(@(~, y) circadian_drift(y, params), t, y0, odeOpts);
X = Y;
end

function [t, X, ouGridSize] = simulate_stochastic_with_ode45(params, y0, totalTime, nSteps, sigmaVec, opts)
noiseClass = normalize_noise_class(opts.noiseClass, 3);
if ~all(strcmpi(noiseClass, 'o'))
    error('Circadian_SDE:UnsupportedNoiseClass', ...
        'The ou_ode45 solver only supports OU noise classes.');
end

t = linspace(0, totalTime, nSteps + 1).';
tau = 1;
xiGrid = utils.sde.generate_ou_path(t, sigmaVec, tau);
odeOpts = odeset('RelTol', opts.odeRelTol, 'AbsTol', opts.odeAbsTol);
[~, Y] = ode45(@(time, y) circadian_rhs_with_ou_bridge(time, y, params, t, xiGrid, tau), t, y0, odeOpts);
X = Y;
ouGridSize = size(xiGrid);
end

function [t, X, ouGridSize] = simulate_stochastic_with_milstein(params, y0, totalTime, nSteps, sigmaVec, opts)
noiseClass = normalize_noise_class(opts.noiseClass, 3);
diffusion = { ...
    @(~, y) diffusion_m(y, params), ...
    @(~, y) diffusion_pc(y, params), ...
    @(~, y) diffusion_pn(y, params)};
diffusionJac = { ...
    @(~, y) diffusion_m_jac(y, params), ...
    @(~, y) diffusion_pc_jac(y, params), ...
    @(~, y) diffusion_pn_jac(y, params)};

[t, Xraw] = utils.sde.milstein(@(~, y) circadian_drift(y, params), diffusion, diffusionJac, ...
    totalTime, nSteps, y0, sigmaVec, noiseClass);
X = Xraw.';
ouGridSize = [nSteps + 1, numel(sigmaVec)];
end

function out = circadian_rhs_with_ou_bridge(time, y, params, tGrid, xiGrid, tau)
xi = utils.sde.evaluate_ou_bridge_mean(time, tGrid, xiGrid, tau);
out = circadian_drift_with_ou(y, params, xi);
end

function out = circadian_drift(y, params)
[Af, ~] = free_active(y(3), params);
beta = 0.1572;

out = zeros(3, 1);
out(1) = beta * (Af / params(2) - y(1));
out(2) = beta * (y(1) - y(2));
out(3) = beta * (y(2) - y(3));
end

function out = circadian_drift_with_ou(y, params, xi)
[Af, ~] = free_active(y(3), params);
beta = 0.1572;
xi = reshape(xi, [], 1);

out = zeros(3, 1);
out(1) = beta * (1 + xi(1)) * (Af / params(2) - y(1));
out(2) = beta * (1 + xi(2)) * (y(1) - y(2));
out(3) = beta * (1 + xi(3)) * (y(2) - y(3));
end

function out = diffusion_m(y, params)
[Af, ~] = free_active(y(3), params);
beta = 0.1572;

out = zeros(3, 1);
out(1) = beta * (Af / params(2) - y(1));
end

function out = diffusion_m_jac(y, params)
[~, dAfdPn] = free_active(y(3), params);
beta = 0.1572;

out = zeros(3, 3);
out(1, 1) = -beta;
out(1, 3) = beta * dAfdPn / params(2);
end

function out = diffusion_pc(y, ~)
beta = 0.1572;

out = zeros(3, 1);
out(2) = beta * (y(1) - y(2));
end

function out = diffusion_pc_jac(~, ~)
beta = 0.1572;

out = zeros(3, 3);
out(2, 1) = beta;
out(2, 2) = -beta;
end

function out = diffusion_pn(y, ~)
beta = 0.1572;

out = zeros(3, 1);
out(3) = beta * (y(2) - y(3));
end

function out = diffusion_pn_jac(~, ~)
beta = 0.1572;

out = zeros(3, 3);
out(3, 2) = beta;
out(3, 3) = -beta;
end

function sigmaVec = normalize_sigma(sigma, expectedNum)
sigmaVec = reshape(double(sigma), 1, []);
if isscalar(sigmaVec)
    sigmaVec = repmat(sigmaVec, 1, expectedNum);
elseif numel(sigmaVec) ~= expectedNum
    error('Circadian_SDE:InvalidSigma', ...
        'Expected sigma to be scalar or length %d.', expectedNum);
end
end

function noiseClass = normalize_noise_class(noiseClassInput, expectedNum)
if ischar(noiseClassInput) || isstring(noiseClassInput)
    noiseClass = cellstr(noiseClassInput);
else
    noiseClass = noiseClassInput;
end

if isscalar(noiseClass)
    noiseClass = repmat(noiseClass, 1, expectedNum);
end
if numel(noiseClass) ~= expectedNum
    error('Circadian_SDE:InvalidNoiseClass', ...
        'Expected noiseClass to be scalar or length %d.', expectedNum);
end
noiseClass = cellfun(@char, noiseClass, 'UniformOutput', false);
end

function [value, derivativePn] = free_active(Pn, params)
Kd = params(1);
AT = params(2);
delta = AT - Pn - Kd;
rootTerm = sqrt(delta.^2 + 4 * Kd * AT);
value = 0.5 * (delta + rootTerm);
derivativePn = 0.5 * (-1 - delta ./ rootTerm);
end

function ensure_shared_sde_path()
persistent isReady
if ~isempty(isReady) && isReady
    return
end

sdeDir = fileparts(mfilename('fullpath'));
circadianDir = fileparts(sdeDir);
fmamDir = fileparts(circadianDir);
addpath(fmamDir);

isReady = true;
end
