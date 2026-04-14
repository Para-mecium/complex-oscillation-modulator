function data = generate_extended_fig2d_data()
scriptDir = fileparts(mfilename('fullpath'));
cacheFile = fullfile(scriptDir, 'extended_fig2d_sde_data.mat');
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
data.meta.nReplicates = 100;
data.meta.distributionSeeds = struct();
data.meta.distributionSeeds.before = 1001:1100;
data.meta.distributionSeeds.after = 1001:1100;
data.meta.distributionOptions = struct();
data.meta.distributionOptions.warmupT = 600;
data.meta.distributionOptions.measureT = 60;
data.meta.distributionOptions.dt = 0.01;
data.meta.distributionOptions.sigma = 0.1;
data.meta.distributionOptions.noiseClass = repmat({'o'}, 1, 8);

data.before = struct();
data.before.distribution = build_distribution(beforeParams, beforeY0, ...
    data.meta.distributionSeeds.before, data.meta.distributionOptions);
data.after = struct();
data.after.distribution = build_distribution(afterParams, afterY0, ...
    data.meta.distributionSeeds.after, data.meta.distributionOptions);

save(cacheFile, '-struct', 'data');
end

function distribution = build_distribution(params, y0, seeds, opts)
runOpts = opts;
runOpts.T = opts.warmupT + opts.measureT;
minS = zeros(numel(seeds), 1);
minH = zeros(numel(seeds), 1);

for i = 1:numel(seeds)
    disp(['Simulating ' num2str(i) 'th seed...'])
    runOpts.seed = seeds(i);
    sim = Longevity_SDE(params, y0, runOpts);
    XMeasure = sim.X(sim.t >= opts.warmupT, :);
    minS(i) = min(XMeasure(:, 3));
    minH(i) = min(XMeasure(:, 4));
end

distribution = struct('minS', minS, 'minH', minH);
end
