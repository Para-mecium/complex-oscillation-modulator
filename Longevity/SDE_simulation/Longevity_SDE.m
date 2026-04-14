function result = Longevity_SDE(params, y0, opts)
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
    opts.noiseClass = repmat({'o'}, 1, 8);
end

params = reshape(double(params), 1, []);
y0 = reshape(double(y0), [], 1);
N = round(opts.T / opts.dt);
rng(double(opts.seed), 'twister');

drift = @(t, y) longevity_drift(y, params);
diffusion = { ...
    @(t, y) diffusion_alphaS(y, params), ...
    @(t, y) diffusion_deltamS(y, params), ...
    @(t, y) diffusion_alphaH(y, params), ...
    @(t, y) diffusion_deltamH(y, params), ...
    @(t, y) diffusion_betaS(y, params), ...
    @(t, y) diffusion_deltaS(y, params), ...
    @(t, y) diffusion_betaH(y, params), ...
    @(t, y) diffusion_deltaH(y, params)};
diffusionJac = { ...
    @(t, y) diffusion_alphaS_jac(y, params), ...
    @(t, y) diffusion_deltamS_jac(y, params), ...
    @(t, y) diffusion_alphaH_jac(y, params), ...
    @(t, y) diffusion_deltamH_jac(y, params), ...
    @(t, y) diffusion_betaS_jac(y, params), ...
    @(t, y) diffusion_deltaS_jac(y, params), ...
    @(t, y) diffusion_betaH_jac(y, params), ...
    @(t, y) diffusion_deltaH_jac(y, params)};

[t, Xraw] = utils.sde.milstein(drift, diffusion, diffusionJac, opts.T, N, y0, opts.sigma, opts.noiseClass);
X = Xraw.';
t = t(:);
if size(X, 1) ~= numel(t) || size(X, 2) ~= 4
    error('Longevity_SDE:InvalidOutputShape', ...
        'Expected X to have size [numel(t), 4].');
end

result = struct();
result.t = t;
result.X = X;
result.params = params;
result.y0 = y0;
result.opts = opts;
result.meta = struct('seed', double(opts.seed));
end

function ensure_shared_sde_path()
persistent isReady
if ~isempty(isReady) && isReady
    return
end

sdeDir = fileparts(mfilename('fullpath'));
longevityDir = fileparts(sdeDir);
fmamDir = fileparts(longevityDir);
addpath(fmamDir);

isReady = true;
end

function out = longevity_drift(y, params)
out = zeros(4, 1);
mS = y(1);
mH = y(2);
S = y(3);
H = y(4);

out(1) = params(3) + params(1) * H.^params(12) ./ (params(10)^params(12) + H.^params(12)) - params(7) * mS;
out(2) = params(4) + params(2) * params(11)^params(13) ./ (params(11)^params(13) + S.^params(13)) - params(7) * mH;
out(3) = params(5) * mS - params(8) * S;
out(4) = params(6) * mH - params(9) * H;
end

function out = diffusion_alphaS(y, params)
out = zeros(4, 1);
H = y(4);
out(1) = params(1) * H.^params(12) ./ (params(10)^params(12) + H.^params(12));
end

function out = diffusion_alphaS_jac(y, params)
out = zeros(4, 4);
H = y(4);
out(1, 4) = params(1) * params(12) * H.^(params(12) - 1) ./ ...
    (params(10)^params(12) + H.^params(12)).^2;
end

function out = diffusion_deltamS(y, params)
out = zeros(4, 1);
out(1) = -params(7) * y(1);
end

function out = diffusion_deltamS_jac(~, params)
out = zeros(4, 4);
out(1, 1) = -params(7);
end

function out = diffusion_alphaH(y, params)
out = zeros(4, 1);
S = y(3);
out(2) = params(2) * params(11)^params(13) ./ (params(11)^params(13) + S.^params(13));
end

function out = diffusion_alphaH_jac(y, params)
out = zeros(4, 4);
S = y(3);
out(2, 3) = params(2) * params(13) * params(11)^(params(13) - 1) ./ ...
    (params(11)^params(13) + S.^params(13)).^2;
end

function out = diffusion_deltamH(y, params)
out = zeros(4, 1);
out(2) = -params(7) * y(2);
end

function out = diffusion_deltamH_jac(~, params)
out = zeros(4, 4);
out(2, 2) = -params(7);
end

function out = diffusion_betaS(y, params)
out = zeros(4, 1);
out(3) = params(5) * y(1);
end

function out = diffusion_betaS_jac(~, params)
out = zeros(4, 4);
out(3, 1) = params(5);
end

function out = diffusion_deltaS(y, params)
out = zeros(4, 1);
out(3) = -params(8) * y(3);
end

function out = diffusion_deltaS_jac(~, params)
out = zeros(4, 4);
out(3, 3) = -params(8);
end

function out = diffusion_betaH(y, params)
out = zeros(4, 1);
out(4) = params(6) * y(2);
end

function out = diffusion_betaH_jac(~, params)
out = zeros(4, 4);
out(4, 2) = params(6);
end

function out = diffusion_deltaH(y, params)
out = zeros(4, 1);
out(4) = -params(9) * y(4);
end

function out = diffusion_deltaH_jac(~, params)
out = zeros(4, 4);
out(4, 4) = -params(9);
end
