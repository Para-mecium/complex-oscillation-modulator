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
    opts.sigma = 0.1;
end
if ~isfield(opts, 'seed')
    opts.seed = 1;
end
if ~isfield(opts, 'noiseClass')
    opts.noiseClass = repmat({'o'}, 1, 3);
end

params = reshape(double(params), 1, []);
y0 = reshape(double(y0), [], 1);
if numel(params) ~= 2
    error('Circadian_SDE:InvalidParams', ...
        'Expected params to be [Kd, AT].');
end
if numel(y0) ~= 3
    error('Circadian_SDE:InvalidInitialState', ...
        'Expected y0 to contain [M; Pc; Pn].');
end

N = round(opts.T / opts.dt);
rng(double(opts.seed), 'twister');

drift = @(t, y) circadian_drift(y, params);
diffusion = { ...
    @(t, y) diffusion_m(y, params), ...
    @(t, y) diffusion_pc(y, params), ...
    @(t, y) diffusion_pn(y, params)};
diffusionJac = { ...
    @(t, y) diffusion_m_jac(y, params), ...
    @(t, y) diffusion_pc_jac(y, params), ...
    @(t, y) diffusion_pn_jac(y, params)};

[t, Xraw] = utils.sde.milstein(drift, diffusion, diffusionJac, opts.T, N, y0, opts.sigma, opts.noiseClass);
X = Xraw.';
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
result.meta = struct('seed', double(opts.seed));
end

function out = circadian_drift(y, params)
[Af, ~] = free_active(y(3), params);
beta = 0.1572;

out = zeros(3, 1);
out(1) = beta * (Af / params(2) - y(1));
out(2) = beta * (y(1) - y(2));
out(3) = beta * (y(2) - y(3));
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
refactorDir = fileparts(sdeDir);
fmamDir = fileparts(refactorDir);
addpath(fmamDir);

isReady = true;
end
