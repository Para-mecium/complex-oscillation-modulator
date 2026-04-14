function tests = test_circadian_fig5d_sde
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(testDir);
addpath(rootDir, '-begin');
addpath(fullfile(rootDir, 'Circadian'), '-begin');
testCase.TestData.rootDir = rootDir;
end

function setup(testCase)
circadian.ensure_paths();

sdeDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Circadian', 'SDE_simulation');
addpath(sdeDir);

stubDir = tempname;
mkdir(stubDir);
write_run_modulation_task_stub(stubDir);
addpath(stubDir, '-begin');
clear run_modulation_task

testCase.TestData.stubDir = stubDir;
testCase.TestData.sdeDir = sdeDir;
testCase.TestData.cleanup = onCleanup(@() cleanup_paths(stubDir, sdeDir));
end

function teardown(~)
close all force
end

function testCircadianSDEUsesStableSeedPolicy(testCase)
params = [8.5e-5, 0.045];
y0 = [0.10; 0.08; 0.07];
opts = struct('T', 4, 'dt', 0.1, 'seed', 11);

first = Circadian_SDE(params, y0, opts);
second = Circadian_SDE(params, y0, opts);
third = Circadian_SDE(params, y0, struct('T', 4, 'dt', 0.1, 'seed', 12));

verifyEqual(testCase, first.X, second.X, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, norm(first.X - third.X, 'fro'), 1e-8);
verifySize(testCase, first.X, [numel(first.t), 3]);
verifyEqual(testCase, first.meta.seed, 11);
verifyEqual(testCase, first.opts.sigma, 0.4, 'AbsTol', 1e-12);
verifyEqual(testCase, first.opts.noiseClass, repmat({'o'}, 1, 3));
verifyEqual(testCase, first.meta.solver, 'ou_ode45');
end

function testCircadianSDEZeroSigmaUsesDeterministicOde45(testCase)
params = [8.5e-5, 0.045];
y0 = [0.10; 0.08; 0.07];

first = Circadian_SDE(params, y0, struct('T', 4, 'dt', 0.1, 'sigma', 0, 'seed', 11));
second = Circadian_SDE(params, y0, struct('T', 4, 'dt', 0.1, 'sigma', 0, 'seed', 12));

verifyEqual(testCase, first.X, second.X, 'AbsTol', 1e-12);
verifyEqual(testCase, first.t, second.t, 'AbsTol', 1e-12);
verifySize(testCase, first.X, [numel(first.t), 3]);
verifyEqual(testCase, first.meta.solver, 'ode45_deterministic');
end

function testCircadianSDESupportsLegacyMilsteinSolver(testCase)
params = [8.5e-5, 0.045];
y0 = [0.10; 0.08; 0.07];
baseOpts = struct('T', 4, 'dt', 0.1, 'sigma', 0.4, 'seed', 11);

modern = Circadian_SDE(params, y0, baseOpts);
legacyOpts = baseOpts;
legacyOpts.solver = 'milstein_ou';
legacy = Circadian_SDE(params, y0, legacyOpts);

verifyEqual(testCase, size(modern.X), size(legacy.X));
verifyEqual(testCase, modern.t, legacy.t, 'AbsTol', 1e-12);
verifyEqual(testCase, legacy.meta.solver, 'milstein_ou');
verifyEqual(testCase, modern.meta.ouGridSize, legacy.meta.ouGridSize);
end

function testOUHelpersAreDeterministicAndShapeStable(testCase)
tGrid = linspace(0, 4, 41).';
sigma = [0.4, 0.4, 0.4];

rng(7, 'twister');
first = utils.sde.generate_ou_path(tGrid, sigma, 1);
rng(7, 'twister');
second = utils.sde.generate_ou_path(tGrid, sigma, 1);
zeroPath = utils.sde.generate_ou_path(tGrid, zeros(1, 3), 1);
bridgeValue = utils.sde.evaluate_ou_bridge_mean(1.25, tGrid, first, 1);

verifyEqual(testCase, first, second, 'AbsTol', 1e-12);
verifyEqual(testCase, zeroPath, zeros(size(zeroPath)), 'AbsTol', 1e-12);
verifySize(testCase, first, [numel(tGrid), 3]);
verifySize(testCase, bridgeValue, [3, 1]);
end

function testOuOde45PeakIsCloserThanLegacyMilstein(testCase)
mark = make_synthetic_mark_data();
mark = mark.markResults{2};
y0 = reshape(mark.orbit.y(1, :), [], 1);
baseOpts = struct( ...
    'T', 1200, ...
    'dt', 0.1, ...
    'sigma', 0.4, ...
    'seed', 13, ...
    'noiseClass', {{'o', 'o', 'o'}}, ...
    'odeRelTol', 1e-8, ...
    'odeAbsTol', 1e-10);

modern = Circadian_SDE(mark.params, y0, baseOpts);
legacyOpts = baseOpts;
legacyOpts.solver = 'milstein_ou';
legacy = Circadian_SDE(mark.params, y0, legacyOpts);
refScaled = 24 / mark.measures.period;
modernPeak = compute_scaled_peak(modern.t, modern.X(:, 2) + modern.X(:, 3), baseOpts.dt);
legacyPeak = compute_scaled_peak(legacy.t, legacy.X(:, 2) + legacy.X(:, 3), baseOpts.dt);

verifyLessThanOrEqual(testCase, abs(modernPeak - refScaled), abs(legacyPeak - refScaled) + 1e-12);
end

function testBuildFig5dSdeDataProducesExpectedStatistics(testCase)
cfg = make_small_cfg();
markData = make_synthetic_mark_data();

data = build_fig5d_sde_data(cfg, markData);
expectedCacheFile = fullfile(cfg.io.cacheDir, 'fig5d', 'sde', ...
    sprintf('data_sigma_%s.mat', circadian.format_sigma_tag(cfg.fig5d.sde.plotSigma)));

verifyEqual(testCase, numel(data.representative), 3);
verifyEqual(testCase, numel(data.psd), 3);
verifyEqual(testCase, numel(data.distribution), 3);
verifyEqual(testCase, data.meta.nReplicates, 3);
verifyEqual(testCase, data.meta.psdSeeds, data.meta.distributionSeeds);
verifyEqual(testCase, numel(data.distribution(1).maxPtot), cfg.fig5d.sde.nReplicates);
verifyEqual(testCase, data.representative(1).seed, cfg.fig5d.sde.representativeSeeds(1));
verifyEqual(testCase, data.meta.solver, 'ou_ode45');
verifyEqual(testCase, data.meta.sigma, cfg.fig5d.sde.plotSigma, 'AbsTol', 1e-12);
verifyEqual(testCase, data.meta.plotSigma, cfg.fig5d.sde.plotSigma, 'AbsTol', 1e-12);
verifyEqual(testCase, data.meta.parallelWorkersRequested, cfg.fig5d.sde.parallelWorkers);
verifyTrue(testCase, isnumeric(data.meta.parallelWorkersUsed) && isscalar(data.meta.parallelWorkersUsed));
verifyTrue(testCase, isfield(data.psd(1), 'frequency'));
verifyTrue(testCase, isfield(data.psd(1), 'mean'));
verifyTrue(testCase, isfield(data.psd(1), 'lower'));
verifyTrue(testCase, isfield(data.psd(1), 'upper'));
verifyTrue(testCase, isfile(expectedCacheFile));
end

function testBuildFig5dSdeDataKeepsSigmaSpecificCachesSeparate(testCase)
cfg = make_small_cfg();
markData = make_synthetic_mark_data();
cfgAlt = cfg;
cfgAlt.fig5d.sde.plotSigma = 0.05;
cfgAlt.fig5d.sde.sigma = 0.05;

build_fig5d_sde_data(cfg, markData);
build_fig5d_sde_data(cfgAlt, markData);

cacheFileA = fullfile(cfg.io.cacheDir, 'fig5d', 'sde', ...
    sprintf('data_sigma_%s.mat', circadian.format_sigma_tag(cfg.fig5d.sde.plotSigma)));
cacheFileB = fullfile(cfg.io.cacheDir, 'fig5d', 'sde', ...
    sprintf('data_sigma_%s.mat', circadian.format_sigma_tag(cfgAlt.fig5d.sde.plotSigma)));

verifyNotEqual(testCase, cacheFileA, cacheFileB);
verifyTrue(testCase, isfile(cacheFileA));
verifyTrue(testCase, isfile(cacheFileB));
end

function testReproduceFig5dBuildsFourPanelsAndSdeBundle(testCase)
cfg = make_small_cfg();

result = reproduce_fig5d(cfg);
titles = get_panel_titles(result.figure);
figureBase = fullfile(cfg.io.figureDir, ...
    sprintf('fig5d_sigma_%s', circadian.format_sigma_tag(cfg.fig5d.sde.plotSigma)));

verifyTrue(testCase, isfield(result.data, 'sde'));
verifyEqual(testCase, numel(titles), 4);
verifyTrue(testCase, any(strcmp(titles, 'Fig. 5d1: Iso-maximum curves')));
verifyTrue(testCase, any(strcmp(titles, 'Fig. 5d2: Representative noisy time series')));
verifyTrue(testCase, any(startsWith(titles, 'Fig. 5d3: Mean PSD')));
verifyTrue(testCase, any(startsWith(titles, 'Fig. 5d4: Maximum distribution')));
verifyTrue(testCase, isfile(fullfile(cfg.io.cacheDir, 'fig5d', 'sde', ...
    sprintf('data_sigma_%s.mat', circadian.format_sigma_tag(cfg.fig5d.sde.plotSigma)))));
verifyEqual(testCase, result.data.figureBase, figureBase);
verifyTrue(testCase, isfile([figureBase '.png']));
verifyTrue(testCase, isfile([figureBase '.pdf']));
end

function titles = get_panel_titles(fig)
axesHandles = findall(fig, 'Type', 'axes');
titles = {};
for i = 1:numel(axesHandles)
    titleText = get(get(axesHandles(i), 'Title'), 'String');
    if ischar(titleText) && startsWith(titleText, 'Fig. 5d')
        titles{end + 1} = titleText; %#ok<AGROW>
    end
end
end

function peakScaled = compute_scaled_peak(~, samples, dt)
signal = reshape(samples, [], 1) - mean(samples);
n = numel(signal);
fs = 1 / dt;
fftValues = fft(signal);
powerValues = abs(fftValues).^2 / (fs * max(n, 1));
halfIdx = floor(n / 2) + 1;
frequency = (0:halfIdx - 1).' * (fs / n);
spectrum = powerValues(1:halfIdx);
if numel(spectrum) > 2
    spectrum(2:end-1) = 2 * spectrum(2:end-1);
end
scaledFrequency = 24 * frequency;
visibleMask = scaledFrequency >= 0.9 & scaledFrequency <= 1.1;
visibleIdx = find(visibleMask);
if isempty(visibleIdx)
    error('test:NoVisiblePsdPeak', ...
        'No PSD samples fall inside the visible frequency band.');
end
[~, peakRelIdx] = max(spectrum(visibleIdx));
peakScaled = scaledFrequency(visibleIdx(peakRelIdx));
end

function cfg = make_small_cfg()
cfg = default_config();
cfg.fig5d.maxima = 0.12;
cfg.fig5d.ATseed = 0.045;
cfg.fig5d.Kdseed = 5e-5;
cfg.fig5d.markMaximum = 0.12;
cfg.fig5d.markPeriods = [23.5, 24, 24.5];
cfg.fig5d.markSeeds = [9e-5, 0.05; 8e-5, 0.045; 7e-5, 0.04];
cfg.fig5d.KdRange = [1e-5, 1.4e-4];
cfg.fig5d.KdAxis = [1e-5, 1.4e-4];
cfg.fig5d.ATAxis = [0.01, 0.09];
cfg.fig5d.maxParamJump = [1, 1];
cfg.fig5d.maxPeriodJump = 10;
cfg.fig5d.sde.representativeSeeds = [1, 2, 3];
cfg.fig5d.sde.psdSeeds = 101:103;
cfg.fig5d.sde.distributionSeeds = 101:103;
cfg.fig5d.sde.nReplicates = 3;
cfg.fig5d.sde.representativeWarmupT = 2;
cfg.fig5d.sde.representativeMeasureT = 4;
cfg.fig5d.sde.psdWarmupT = 2;
cfg.fig5d.sde.psdMeasureT = 4;
cfg.fig5d.sde.distributionWarmupT = 2;
cfg.fig5d.sde.distributionMeasureT = 4;
cfg.fig5d.sde.dt = 0.1;
cfg.fig5d.sde.plotSigma = 0.4;
cfg.fig5d.sde.sigma = 0.4;
cfg.fig5d.sde.parallelWorkers = 1;
cfg.fig5d.sde.solver = 'ou_ode45';
cfg.fig5d.sde.odeRelTol = 1e-8;
cfg.fig5d.sde.odeAbsTol = 1e-10;
cfg.fig5d.sde.psdBandPercentiles = [10, 90];

tempRoot = tempname;
mkdir(tempRoot);
cfg.io.cacheDir = fullfile(tempRoot, 'cache');
cfg.io.figureDir = fullfile(tempRoot, 'figure');
end

function markData = make_synthetic_mark_data()
periods = [23.5, 24, 24.5];
params = [1.4e-4, 0.051; 8.34e-5, 0.046; 4.43e-5, 0.043];
marks = cell(1, numel(periods));

for i = 1:numel(periods)
    t = linspace(0, periods(i), 25).';
    ptot = 0.12 - 0.02 + 0.02 * (1 + sin(2 * pi * t / periods(i)));
    pc = 0.58 * ptot;
    pn = 0.42 * ptot;
    marks{i} = struct( ...
        'params', params(i, :), ...
        'measures', struct('period', periods(i), 'maximum', 0.12), ...
        'orbit', struct('t', t, 'y', [0.8 * ptot, pc, pn], 'obs', ptot, 'period', periods(i)));
end

markData = struct('markResults', {marks}, 'markPeriods', periods);
end

function cleanup_paths(stubDir, sdeDir)
if isfolder(stubDir)
    rmpath(stubDir);
    clear run_modulation_task
    try
        rmdir(stubDir, 's');
    catch
    end
end

if isfolder(sdeDir)
    rmpath(sdeDir);
end
end

function write_run_modulation_task_stub(stubDir)
stubPath = fullfile(stubDir, 'run_modulation_task.m');
fid = fopen(stubPath, 'w');
assert(fid >= 0, 'Failed to create run_modulation_task stub.');
cleanup = onCleanup(@() fclose(fid));

lines = { ...
    'function result = run_modulation_task(cfg)', ...
    'goalOrder = cfg.goalOrder;', ...
    'if isequal(goalOrder, {''maximum'', ''Kd''})', ...
    '    result = make_maximum_kd_branch(cfg);', ...
    'elseif isequal(goalOrder, {''maximum'', ''period''})', ...
    '    params = reshape(cfg.startParams, 1, []);', ...
    '    amplitude = max(0.5 * cfg.goals.maximum, 0.01);', ...
    '    result = make_result({''Kd'', ''AT''}, params, cfg.goals.period, amplitude, cfg.goals.maximum, []);', ...
    'else', ...
    '    error(''test:UnexpectedGoalOrder'', ''Unexpected goalOrder in stub.'');', ...
    'end', ...
    'end', ...
    '', ...
    'function result = make_maximum_kd_branch(cfg)', ...
    'params = resolve_start_params(cfg);', ...
    'targetKd = cfg.goals.Kd;', ...
    'atValues = [params(2); params(2) + 0.01 * sign(targetKd - params(1) + eps)];', ...
    'period = [cfg.goals.maximum * 100 + 12; cfg.goals.maximum * 100 + 13];', ...
    'amplitude = max(0.5 * cfg.goals.maximum, 0.01) * [1; 1.05];', ...
    'maximum = cfg.goals.maximum * [1; 1];', ...
    'stopReason = ''target_reached'';', ...
    'paramMatrix = [params(1), atValues(1); targetKd, atValues(2)];', ...
    'result = make_result({''Kd'', ''AT''}, paramMatrix(end, :), period(end), amplitude(end), maximum(end), make_path({''Kd'', ''AT''}, paramMatrix, period, amplitude, maximum, stopReason));', ...
    'end', ...
    '', ...
    'function params = resolve_start_params(cfg)', ...
    'if isfield(cfg, ''startParams'') && ~isempty(cfg.startParams)', ...
    '    params = reshape(cfg.startParams, 1, []);', ...
    'elseif isfield(cfg, ''initialSolverView'') && isfield(cfg.initialSolverView, ''params'')', ...
    '    params = reshape(cfg.initialSolverView.params, 1, []);', ...
    'else', ...
    '    error(''test:MissingStartParams'', ''Missing start params in stub.'');', ...
    'end', ...
    'end', ...
    '', ...
    'function result = make_result(paramNames, params, period, amplitude, maximum, path)', ...
    'orbit = make_orbit(period, params, amplitude, maximum);', ...
    'result = struct();', ...
    'result.params = reshape(params, 1, []);', ...
    'result.orbit = orbit;', ...
    'result.path = path;', ...
    'result.measures = struct(''period'', period, ''amplitude'', amplitude, ''maximum'', maximum);', ...
    'result.solverView = struct(''params'', reshape(params, 1, []));', ...
    'result.derivedView = struct(''period'', period, ''obsAmp'', amplitude, ''obsMax'', maximum, ''obsMin'', maximum - 2 * amplitude);', ...
    'end', ...
    '', ...
    'function path = make_path(paramNames, paramMatrix, period, amplitude, maximum, stopReason)', ...
    'path = struct();', ...
    'path.paramNames = paramNames;', ...
    'path.paramMatrix = paramMatrix;', ...
    'path.period = period(:);', ...
    'path.obsAmplitude = amplitude(:);', ...
    'path.obsMax = maximum(:);', ...
    'path.obsMin = path.obsMax - 2 * path.obsAmplitude;', ...
    'path.lambda = linspace(0, 1, size(paramMatrix, 1)).'';', ...
    'path.directConditionEstimate = ones(size(paramMatrix, 1), 1);', ...
    'path.finalConditionEstimate = ones(size(paramMatrix, 1), 1);', ...
    'path.success = true(size(paramMatrix, 1), 1);', ...
    'path.stopReason = stopReason;', ...
    'path.stopLambda = 1;', ...
    'path.stopTriggerValue = NaN;', ...
    'end', ...
    '', ...
    'function orbit = make_orbit(period, params, amplitude, maximum)', ...
    't = linspace(0, period, 8).'';', ...
    'ptot = maximum - amplitude + amplitude * (1 + sin(2 * pi * t / max(period, eps)));', ...
    'pc = 0.6 * ptot;', ...
    'pn = 0.4 * ptot;', ...
    'orbit = struct();', ...
    'orbit.t = t;', ...
    'orbit.y = [0.8 * ptot, pc, pn];', ...
    'orbit.obs = ptot;', ...
    'orbit.period = period;', ...
    'orbit.params = params;', ...
    'orbit.amplitude = amplitude;', ...
    'orbit.yMax = maximum;', ...
    'orbit.yMin = maximum - 2 * amplitude;', ...
    'end'};
fprintf(fid, '%s\n', lines{:});
end
