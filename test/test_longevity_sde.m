function tests = test_longevity_sde
tests = functiontests(localfunctions);
end

function setup(testCase)
sdeDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Longevity', 'SDE_simulation');
addpath(sdeDir);
testCase.TestData.sdeDir = sdeDir;
testCase.TestData.cleanup = onCleanup(@() rmpath_safe(sdeDir));
end

function testLongevitySDEUsesStableSeedPolicy(testCase)
params = baseline_params();
y0 = baseline_y0();
opts = struct('T', 4, 'dt', 0.1, 'seed', 11);

first = Longevity_SDE(params, y0, opts);
second = Longevity_SDE(params, y0, opts);
third = Longevity_SDE(params, y0, struct('T', 4, 'dt', 0.1, 'seed', 12));

verifyEqual(testCase, first.X, second.X, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, norm(first.X - third.X, 'fro'), 1e-8);
verifySize(testCase, first.X, [numel(first.t), 4]);
verifyEqual(testCase, first.meta.seed, 11);
verifyEqual(testCase, first.opts.sigma, 0.1, 'AbsTol', 1e-12);
verifyEqual(testCase, first.opts.noiseClass, repmat({'o'}, 1, 8));
end

function testLongevitySDEAcceptsVectorSigmaAndScalarNoiseClass(testCase)
params = baseline_params();
y0 = baseline_y0();
opts = struct( ...
    'T', 4, ...
    'dt', 0.1, ...
    'seed', 7, ...
    'sigma', 0.01 * (1:8), ...
    'noiseClass', 'o');

result = Longevity_SDE(params, y0, opts);

verifySize(testCase, result.X, [numel(result.t), 4]);
verifyEqual(testCase, result.opts.sigma, opts.sigma, 'AbsTol', 1e-12);
verifyEqual(testCase, result.opts.noiseClass, 'o');
end

function params = baseline_params()
params = [30.5, 183, 0.1, 0.1, 3.7, 3.7, 0.3, 0.2, 3.8, 326, 185, 3, 4.8];
end

function y0 = baseline_y0()
y0 = [23.3896986043691; 222.452142059269; 207.783959243676; 216.598138352976];
end

function rmpath_safe(pathValue)
if contains(path, pathValue)
    rmpath(pathValue);
end
end
