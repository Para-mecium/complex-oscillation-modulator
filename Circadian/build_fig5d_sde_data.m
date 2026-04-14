function data = build_fig5d_sde_data(cfg, markData)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
sdeCfg = normalize_sde_cfg(cfg.fig5d.sde);
cfg.fig5d.sde = sdeCfg;
cacheFile = fullfile(cfg.io.cacheDir, 'fig5d', 'sde', ...
    sprintf('data_sigma_%s.mat', circadian.format_sigma_tag(sdeCfg.plotSigma)));

if nargin < 2 || isempty(markData)
    data = circadian.cache_get_or_create(cacheFile, @() build_cached_data(cfg));
else
    data = circadian.cache_get_or_create(cacheFile, @() build_uncached_data(cfg, markData));
end
end

function data = build_cached_data(cfg)
markData = build_fig5d_marks(cfg);
data = build_uncached_data(cfg, markData);
end

function data = build_uncached_data(cfg, markData)
sdeDir = fullfile(cfg.io.rootDir, 'SDE_simulation');
if exist('Circadian_SDE', 'file') ~= 2
    addpath(sdeDir);
end

sdeCfg = cfg.fig5d.sde;
markResults = markData.markResults;
nMarks = numel(markResults);
parallelState = prepare_parallel_repeats(sdeCfg);

validate_seed_policy(sdeCfg);

representative = repmat(struct( ...
    'targetPeriod', [], ...
    'seed', [], ...
    'params', [], ...
    't', [], ...
    'X', [], ...
    'Ptot', []), 1, nMarks);
psd = repmat(struct( ...
    'targetPeriod', [], ...
    'params', [], ...
    'seeds', [], ...
    'frequency', [], ...
    'mean', [], ...
    'lower', [], ...
    'upper', []), 1, nMarks);
distribution = repmat(struct( ...
    'targetPeriod', [], ...
    'params', [], ...
    'seeds', [], ...
    'maxPtot', []), 1, nMarks);

for i = 1:nMarks
    mark = markResults{i};
    y0 = reshape(mark.orbit.y(1, :), [], 1);
    log_sde_progress('mark', i, nMarks, mark.measures.period);

    representative(i) = build_representative(mark, y0, sdeCfg, i);
    repeatStats = build_repeat_stats(mark, y0, sdeCfg, parallelState);
    psd(i) = repeatStats.psd;
    distribution(i) = repeatStats.distribution;
end

data = struct();
data.meta = struct( ...
    'nReplicates', sdeCfg.nReplicates, ...
    'representativeSeeds', reshape(double(sdeCfg.representativeSeeds), 1, []), ...
    'psdSeeds', reshape(double(sdeCfg.psdSeeds), 1, []), ...
    'distributionSeeds', reshape(double(sdeCfg.distributionSeeds), 1, []), ...
    'solver', resolve_solver_label(sdeCfg), ...
    'repeatKind', 'independent noise realizations', ...
    'sigma', double(sdeCfg.sigma), ...
    'plotSigma', double(sdeCfg.plotSigma), ...
    'tauNoise', 1, ...
    'noiseClass', {sdeCfg.noiseClass}, ...
    'parallelWorkersRequested', double(parallelState.requestedWorkers), ...
    'parallelWorkersUsed', double(parallelState.usedWorkers), ...
    'psdBandPercentiles', reshape(double(sdeCfg.psdBandPercentiles), 1, []), ...
    'formalStatisticalTests', 'none', ...
    'markPeriods', cellfun(@(mark) mark.measures.period, markResults), ...
    'markMaximum', cfg.fig5d.markMaximum);
data.representative = representative;
data.psd = psd;
data.distribution = distribution;
end

function representative = build_representative(mark, y0, sdeCfg, idx)
runOpts = make_run_opts(sdeCfg, sdeCfg.representativeWarmupT + sdeCfg.representativeMeasureT);
runOpts.seed = sdeCfg.representativeSeeds(idx);
fprintf('[circadian.fig5d.sde] representative: T=%.1f, seed=%d\n', ...
    mark.measures.period, runOpts.seed);

sim = Circadian_SDE(mark.params, y0, runOpts);
measureMask = sim.t >= sdeCfg.representativeWarmupT;
t = sim.t(measureMask) - sdeCfg.representativeWarmupT;
X = sim.X(measureMask, :);

representative = struct( ...
    'targetPeriod', mark.measures.period, ...
    'seed', double(runOpts.seed), ...
    'params', reshape(mark.params, 1, []), ...
    't', t, ...
    'X', X, ...
    'Ptot', X(:, 2) + X(:, 3));
end

function repeatStats = build_repeat_stats(mark, y0, sdeCfg, parallelState)
markParams = reshape(mark.params, 1, []);
psdSeeds = reshape(double(sdeCfg.psdSeeds), 1, []);
distributionSeeds = reshape(double(sdeCfg.distributionSeeds), 1, []);
psdTotalTime = sdeCfg.psdWarmupT + sdeCfg.psdMeasureT;
distributionTotalTime = sdeCfg.distributionWarmupT + sdeCfg.distributionMeasureT;
allSeeds = unique([psdSeeds, distributionSeeds], 'stable');
numSeeds = numel(allSeeds);
totalTimes = zeros(1, numSeeds);
repeatKinds = cell(1, numSeeds);
for i = 1:numSeeds
    usesPsd = any(psdSeeds == allSeeds(i));
    usesDistribution = any(distributionSeeds == allSeeds(i));
    totalTimes(i) = max(usesPsd * psdTotalTime, usesDistribution * distributionTotalTime);
    repeatKinds{i} = resolve_repeat_kind(usesPsd, usesDistribution);
end

simBank = collect_sim_bank(markParams, mark.measures.period, y0, sdeCfg, allSeeds, totalTimes, repeatKinds, parallelState);

repeatStats = struct();
repeatStats.psd = build_psd_from_bank(mark, sdeCfg, psdSeeds, simBank);
repeatStats.distribution = build_distribution_from_bank(mark, sdeCfg, distributionSeeds, simBank);
end

function simBank = collect_sim_bank(markParams, targetPeriod, y0, sdeCfg, allSeeds, totalTimes, repeatKinds, parallelState)
numSeeds = numel(allSeeds);
simBank = repmat(struct('seed', [], 'sim', []), 1, numSeeds);
simCells = cell(1, numSeeds);

if parallelState.enabled && numSeeds > 1
    queue = create_repeat_progress_queue(targetPeriod, numSeeds);
    parfor i = 1:numSeeds
        runOpts = make_run_opts(sdeCfg, totalTimes(i));
        runOpts.seed = allSeeds(i);
        if ~isempty(queue)
            send(queue, struct('kind', repeatKinds{i}, 'index', i, 'seed', runOpts.seed));
        end
        simCells{i} = Circadian_SDE(markParams, y0, runOpts);
    end
else
    for i = 1:numSeeds
        runOpts = make_run_opts(sdeCfg, totalTimes(i));
        runOpts.seed = allSeeds(i);
        log_repeat_progress(repeatKinds{i}, targetPeriod, i, numSeeds, runOpts.seed);
        simCells{i} = Circadian_SDE(markParams, y0, runOpts);
    end
end

for i = 1:numSeeds
    simBank(i).seed = allSeeds(i);
    simBank(i).sim = simCells{i};
end
end

function stats = build_psd_from_bank(mark, sdeCfg, seeds, simBank)
firstPn = select_measurement_window(fetch_pn(simBank, seeds(1)), sdeCfg.psdWarmupT);
[frequency, firstPsd] = single_sided_periodogram(firstPn, sdeCfg.dt);
psdValues = zeros(numel(seeds), numel(firstPsd));
psdValues(1, :) = firstPsd(:).';

for i = 2:numel(seeds)
    Pn = select_measurement_window(fetch_pn(simBank, seeds(i)), sdeCfg.psdWarmupT);
    [~, currentPsd] = single_sided_periodogram(Pn, sdeCfg.dt);
    psdValues(i, :) = currentPsd(:).';
end

percentiles = reshape(double(sdeCfg.psdBandPercentiles), 1, []);
band = prctile(psdValues, percentiles, 1);
stats = struct( ...
    'targetPeriod', mark.measures.period, ...
    'params', reshape(mark.params, 1, []), ...
    'seeds', seeds, ...
    'frequency', frequency(:), ...
    'mean', mean(psdValues, 1).', ...
    'lower', band(1, :).', ...
    'upper', band(2, :).');
end

function stats = build_distribution_from_bank(mark, sdeCfg, seeds, simBank)
maxPtot = zeros(numel(seeds), 1);

for i = 1:numel(seeds)
    Ptot = select_measurement_window(fetch_ptot(simBank, seeds(i)), sdeCfg.distributionWarmupT);
    maxPtot(i) = max(Ptot);
end

stats = struct( ...
    'targetPeriod', mark.measures.period, ...
    'params', reshape(mark.params, 1, []), ...
    'seeds', seeds, ...
    'maxPtot', maxPtot);
end

function signal = fetch_pn(simBank, seed)
matchIdx = find([simBank.seed] == seed, 1);
if isempty(matchIdx)
    error('circadian:MissingSimulationForSeed', ...
        'Missing simulation for seed %d.', seed);
end

sim = simBank(matchIdx).sim;
signal = struct('t', sim.t, 'value', sim.X(:, 3));
end

function signal = fetch_ptot(simBank, seed)
matchIdx = find([simBank.seed] == seed, 1);
if isempty(matchIdx)
    error('circadian:MissingSimulationForSeed', ...
        'Missing simulation for seed %d.', seed);
end

sim = simBank(matchIdx).sim;
signal = struct('t', sim.t, 'value', sim.X(:, 2) + sim.X(:, 3));
end

function values = select_measurement_window(signal, warmupT)
measureMask = signal.t >= warmupT;
values = signal.value(measureMask);
end

function kind = resolve_repeat_kind(usesPsd, usesDistribution)
if usesPsd && usesDistribution
    kind = 'psd+distribution';
elseif usesPsd
    kind = 'psd';
else
    kind = 'distribution';
end
end

function runOpts = make_run_opts(sdeCfg, totalTime)
runOpts = struct();
runOpts.T = totalTime;
runOpts.dt = sdeCfg.dt;
runOpts.sigma = sdeCfg.sigma;
runOpts.noiseClass = sdeCfg.noiseClass;
runOpts.solver = sdeCfg.solver;
runOpts.odeRelTol = sdeCfg.odeRelTol;
runOpts.odeAbsTol = sdeCfg.odeAbsTol;
end

function queue = create_repeat_progress_queue(targetPeriod, total)
queue = [];
if exist('parallel.pool.DataQueue', 'class') ~= 8
    return
end

queue = parallel.pool.DataQueue;
afterEach(queue, @(payload) log_repeat_progress(payload.kind, targetPeriod, total, payload));
end

function state = prepare_parallel_repeats(sdeCfg)
state = struct('enabled', false, 'requestedWorkers', 0, 'usedWorkers', 0);
requestedWorkers = normalize_parallel_workers(sdeCfg.parallelWorkers);
state.requestedWorkers = requestedWorkers;

if requestedWorkers == 1 || ~has_parallel_support()
    log_parallel_setup(requestedWorkers, 0);
    return
end

try
    pool = gcp('nocreate');
    if requestedWorkers > 1
        if isempty(pool) || pool.NumWorkers ~= requestedWorkers
            if ~isempty(pool)
                delete(pool);
            end
            pool = parpool('local', requestedWorkers);
        end
    else
        if isempty(pool)
            pool = parpool('local');
        end
    end

    state.enabled = ~isempty(pool) && pool.NumWorkers > 1;
    if ~isempty(pool)
        state.usedWorkers = pool.NumWorkers;
    end
catch err
    log_parallel_warning(requestedWorkers, err);
    state.enabled = false;
    state.usedWorkers = 0;
end

log_parallel_setup(requestedWorkers, state.usedWorkers);
end

function tf = has_parallel_support()
tf = license('test', 'Distrib_Computing_Toolbox') ...
    && ~isempty(ver('parallel')) ...
    && exist('gcp', 'file') == 2 ...
    && exist('parpool', 'file') == 2;
end

function workers = normalize_parallel_workers(value)
if isempty(value)
    workers = 0;
    return
end

workers = double(value);
if ~(isscalar(workers) && isfinite(workers) && workers >= 0 && workers == round(workers))
    error('circadian:InvalidParallelWorkers', ...
        'cfg.fig5d.sde.parallelWorkers must be empty, zero, or a nonnegative integer.');
end
end

function log_parallel_setup(requestedWorkers, usedWorkers)
if requestedWorkers > 0
    fprintf('[circadian.fig5d.sde] parallel workers: requested=%d, used=%d\n', ...
        requestedWorkers, usedWorkers);
else
    fprintf('[circadian.fig5d.sde] parallel workers: requested=auto, used=%d\n', ...
        usedWorkers);
end
end

function log_parallel_warning(requestedWorkers, err)
if requestedWorkers > 0
    requestedLabel = sprintf('%d', requestedWorkers);
else
    requestedLabel = 'auto';
end

fprintf('[circadian.fig5d.sde] parallel fallback: requested=%s, reason=%s\n', ...
    requestedLabel, err.message);
end

function sdeCfg = normalize_sde_cfg(sdeCfg)
if ~isfield(sdeCfg, 'plotSigma') || isempty(sdeCfg.plotSigma)
    sdeCfg.plotSigma = sdeCfg.sigma;
end
sdeCfg.plotSigma = double(sdeCfg.plotSigma);
sdeCfg.sigma = sdeCfg.plotSigma;

if ~isfield(sdeCfg, 'parallelWorkers')
    sdeCfg.parallelWorkers = [];
end
end

function [frequency, spectrum] = single_sided_periodogram(samples, dt)
samples = reshape(double(samples), [], 1);
samples = samples - mean(samples);

n = numel(samples);
fs = 1 / dt;
fftValues = fft(samples);
powerValues = abs(fftValues).^2 / (fs * max(n, 1));

halfIdx = floor(n / 2) + 1;
frequency = (0:halfIdx - 1).' * (fs / n);
spectrum = powerValues(1:halfIdx);
if numel(spectrum) > 2
    spectrum(2:end-1) = 2 * spectrum(2:end-1);
end
end

function validate_seed_policy(sdeCfg)
if numel(sdeCfg.psdSeeds) ~= sdeCfg.nReplicates
    error('circadian:PsdSeedMismatch', ...
        'cfg.fig5d.sde.psdSeeds must have cfg.fig5d.sde.nReplicates entries.');
end
if numel(sdeCfg.distributionSeeds) ~= sdeCfg.nReplicates
    error('circadian:DistributionSeedMismatch', ...
        'cfg.fig5d.sde.distributionSeeds must have cfg.fig5d.sde.nReplicates entries.');
end
end

function label = resolve_solver_label(sdeCfg)
sigma = double(sdeCfg.sigma);
if all(sigma(:) == 0)
    label = 'ode45_deterministic';
else
    label = char(sdeCfg.solver);
end
end

function log_sde_progress(stage, index, total, targetPeriod)
fprintf('[circadian.fig5d.sde] %s %d/%d: T=%.1f\n', ...
    stage, index, total, targetPeriod);
end

function log_repeat_progress(kind, targetPeriod, indexOrTotal, totalOrPayload, seed)
if nargin == 4 && isstruct(totalOrPayload)
    payload = totalOrPayload;
    index = payload.index;
    total = indexOrTotal;
    seed = payload.seed;
else
    index = indexOrTotal;
    total = totalOrPayload;
end

fprintf('[circadian.fig5d.sde] %s %d/%d: T=%.1f, seed=%d\n', ...
    kind, index, total, targetPeriod, seed);
end
