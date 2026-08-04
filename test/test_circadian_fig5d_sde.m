function tests = test_circadian_fig5d_sde
%TEST_CIRCADIAN_FIG5D_SDE Regression coverage for the current SDE simulator.

tests = functiontests(localfunctions);
end

function setup(testCase)
sdeDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'Circadian', 'SDE_simulation');
addpath(sdeDir);
testCase.TestData.cleanup = onCleanup(@() rmpath_safe(sdeDir));
end

function testCircadianSDEUsesStableSeedPolicy(testCase)
params = baseline_params();
y0 = baseline_y0();
opts = struct('T', 4, 'dt', 0.1, 'seed', 11);

first = Circadian_SDE(params, y0, opts);
second = Circadian_SDE(params, y0, opts);
third = Circadian_SDE(params, y0, struct('T', 4, 'dt', 0.1, 'seed', 12));

verifyEqual(testCase, first.X, second.X, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, norm(first.X - third.X, 'fro'), 1e-8);
verifySize(testCase, first.X, [numel(first.t), 3]);
verifyEqual(testCase, first.meta.seed, 11);
verifyEqual(testCase, first.opts.sigma, 0.1, 'AbsTol', 1e-12);
verifyEqual(testCase, first.opts.noiseClass, repmat({'o'}, 1, 3));
end

function testCircadianSDEZeroSigmaIsSeedIndependent(testCase)
params = baseline_params();
y0 = baseline_y0();

first = Circadian_SDE(params, y0, ...
    struct('T', 4, 'dt', 0.1, 'sigma', 0, 'seed', 11));
second = Circadian_SDE(params, y0, ...
    struct('T', 4, 'dt', 0.1, 'sigma', 0, 'seed', 12));

verifyEqual(testCase, first.X, second.X, 'AbsTol', 1e-12);
verifyEqual(testCase, first.t, second.t, 'AbsTol', 1e-12);
verifySize(testCase, first.X, [numel(first.t), 3]);
end

function testCircadianSDEAcceptsVectorSigmaAndScalarNoiseClass(testCase)
params = baseline_params();
y0 = baseline_y0();
opts = struct( ...
    'T', 4, ...
    'dt', 0.1, ...
    'seed', 7, ...
    'sigma', [0.01, 0.02, 0.03], ...
    'noiseClass', 'o');

result = Circadian_SDE(params, y0, opts);

verifySize(testCase, result.X, [numel(result.t), 3]);
verifyEqual(testCase, result.opts.sigma, opts.sigma, 'AbsTol', 1e-12);
verifyEqual(testCase, result.opts.noiseClass, 'o');
end

function params = baseline_params()
params = [8.5e-5, 0.045];
end

function y0 = baseline_y0()
y0 = [0.10; 0.08; 0.07];
end

function rmpath_safe(pathValue)
if contains(path, pathValue)
    rmpath(pathValue);
end
end
