function tests = test_circuit_forward_orbit
%TEST_CIRCUIT_FORWARD_ORBIT Coverage for the RLT circuit forward-orbit wrapper.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir, '-begin');
addpath(fullfile(rootDir, 'RLT_circuit'), '-begin');
addpath(fullfile(rootDir, 'PO_extract'), '-begin');
testCase.TestData.rootDir = rootDir;
end

function testCircuitForwardOrbitMatchesObservableAndFeatureContracts(testCase)
data = load(fullfile(testCase.TestData.rootDir, 'RLT_circuit', ...
    'parameter_inference', 'initData_ODE.mat'));
Parameters = reshape(data.Parameters, 1, []);
y0 = data.TS{2}(1, :).';

result = circuit_forward_orbit(Parameters, y0, struct());

verifyTrue(testCase, result.success, ...
    sprintf('Expected circuit_forward_orbit to succeed, got: %s', string(result.msg)));
verifyEqual(testCase, result.Parameters, Parameters, 'AbsTol', 1e-12);
verifyEqual(testCase, result.initialState, y0, 'AbsTol', 1e-12);
verifyEqual(testCase, result.features.period, result.orbit.period, 'AbsTol', 1e-10);
verifyEqual(testCase, result.orbit.t(1), 0, 'AbsTol', 1e-12);
verifyEqual(testCase, result.orbit.t(end), result.orbit.period, 'AbsTol', 1e-8);
verifyEqual(testCase, size(result.orbit.y, 2), 3);
verifyGreaterThan(testCase, min(result.features.state.amplitude), 0);
verifyEqual(testCase, result.features.observable.series(:, 1), ...
    result.orbit.y(:, 1) + result.orbit.y(:, 2), 'AbsTol', 1e-12);
verifyLessThanOrEqual(testCase, norm(result.orbit.y(end, :) - result.orbit.y(1, :), inf), 1e-6);
end
