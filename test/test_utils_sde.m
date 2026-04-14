function tests = test_utils_sde
tests = functiontests(localfunctions);
end

function setup(testCase)
fmamDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fmamDir);
testCase.TestData.fmamDir = fmamDir;
testCase.TestData.cleanup = onCleanup(@() rmpath_safe(fmamDir));
end

function testGenerateOuPathAndBridgeAreReproducible(testCase)
tGrid = linspace(0, 4, 41).';
sigma = [0.4, 0.2, 0.1];

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

function testMilsteinWhiteNoiseIsReproducible(testCase)
drift = @(~, x) -0.25 * x;
diffusion = @(~, ~) 1;
diffusionJac = @(~, ~) 0;

rng(9, 'twister');
[t1, x1] = utils.sde.milstein(drift, diffusion, diffusionJac, 2, 40, 0.5, 0.2, 'w');
rng(9, 'twister');
[t2, x2] = utils.sde.milstein(drift, diffusion, diffusionJac, 2, 40, 0.5, 0.2, 'w');

verifyEqual(testCase, t1, t2, 'AbsTol', 1e-12);
verifyEqual(testCase, x1, x2, 'AbsTol', 1e-12);
verifySize(testCase, x1, [1, 41]);
end

function testMilsteinOuNoiseSupportsScalarClassAndVectorSigma(testCase)
drift = @(~, x) -0.1 * x;
diffusion = {@(~, ~) 1, @(~, ~) -2};
diffusionJac = {@(~, ~) 0, @(~, ~) 0};

rng(11, 'twister');
[t1, x1] = utils.sde.milstein(drift, diffusion, diffusionJac, 1, 20, 0.5, [0.2, 0.1], 'o');
rng(11, 'twister');
[t2, x2] = utils.sde.milstein(drift, diffusion, diffusionJac, 1, 20, 0.5, [0.2, 0.1], 'o');

verifyEqual(testCase, t1, t2, 'AbsTol', 1e-12);
verifyEqual(testCase, x1, x2, 'AbsTol', 1e-12);
verifySize(testCase, x1, [1, 21]);
end

function rmpath_safe(pathValue)
if contains(path, pathValue)
    rmpath(pathValue);
end
end
