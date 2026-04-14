function data = generate_extended_fig2c_data()
scriptDir = fileparts(mfilename('fullpath'));
cacheFile = fullfile(scriptDir, 'extended_fig2c_sde_data.mat');
targetCacheFile = fullfile(scriptDir, 'extended_fig2b_cache.mat');

if ~isfile(targetCacheFile)
    generate_extended_fig2b_data();
end

targetCache = load(targetCacheFile);
addpath(fullfile(scriptDir, 'SDE_simulation'));

beforeParams = [30.5, 183, 0.1, 0.1, 3.7, 3.7, 0.3, 0.2, 3.8, 326, 185, 3, 4.8];
beforeY0 = [23.3896986043691; 222.452142059269; 207.783959243676; 216.598138352976];
afterParams = reshape(double(targetCache.target.Parameters), 1, []);
afterY0 = reshape(double(targetCache.target.TS_var(1, :)), [], 1);

data = struct();
data.meta = struct();
data.meta.agingZone = struct('S', 200, 'H', 100);
data.meta.representativeSeeds = struct();
data.meta.representativeSeeds.before = 1;
data.meta.representativeSeeds.after = 1;
data.meta.representativeOptions = struct();
data.meta.representativeOptions.warmupT = 1000;
data.meta.representativeOptions.measureT = 60;
data.meta.representativeOptions.dt = 0.01;
data.meta.representativeOptions.sigma = 0.1;
data.meta.representativeOptions.noiseClass = repmat({'o'}, 1, 8);

repOpts = data.meta.representativeOptions;
repOpts.seed = data.meta.representativeSeeds.before;
before = build_representative(beforeParams, beforeY0, repOpts);
repOpts.seed = data.meta.representativeSeeds.after;
after = build_representative(afterParams, afterY0, repOpts);

data.before = struct();
data.before.representative = struct('t', before.t, 'X', before.X);
data.after = struct();
data.after.representative = struct('t', after.t, 'X', after.X);

save(cacheFile, '-struct', 'data');
end

function representative = build_representative(params, y0, opts)
runOpts = opts;
runOpts.T = opts.warmupT + opts.measureT;
sim = Longevity_SDE(params, y0, runOpts);
measureMask = sim.t >= opts.warmupT;

representative = struct();
representative.t = sim.t(measureMask) - opts.warmupT;
representative.X = sim.X(measureMask, :);
end
