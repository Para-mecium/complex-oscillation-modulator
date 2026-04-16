function tests = test_orbit_feature_wrapper_compatibility
%TEST_ORBIT_FEATURE_WRAPPER_COMPATIBILITY Wrapper compatibility checks.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir, '-begin');
addpath(fullfile(rootDir, 'PO_extract'), '-begin');
addpath(fullfile(rootDir, 'flexible_modulators'), '-begin');
addpath(fullfile(rootDir, 'Circadian'), '-begin');
testCase.TestData.rootDir = rootDir;
end

function testFlexmodOrbitMatchesCoarseFeatures(testCase)
cfg = make_flexmod_cfg();
model = struct('rhs', @harmonic_rhs);

orbit = flexmod.find_orbit(model, 1, cfg, cfg.initialState);

verifyTrue(testCase, orbit.success, 'Expected find_orbit to return a periodic orbit.');
verifyGreaterThan(testCase, orbit.period, 6);
verifyLessThan(testCase, orbit.period, 7);

features = evaluate_orbit_features(struct('t', orbit.t, 'y', orbit.y, 'period', orbit.period), ...
    [], ["period", "state_max", "state_min", "state_amplitude"], ...
    struct('extremaRefinement', false));
verifyEqual(testCase, orbit.amplitude, features.state.amplitude(orbit.proteinIndex), 'AbsTol', 1e-10);
verifyEqual(testCase, orbit.yMax, features.state.max(orbit.proteinIndex), 'AbsTol', 1e-10);
verifyEqual(testCase, orbit.yMin, features.state.min(orbit.proteinIndex), 'AbsTol', 1e-10);
end

function testCircadianOrbitMatchesCoarseObservableFeatures(testCase)
cfg = circadian.default_config();
cfg.initialState = [1; 0; 0];
cfg.orbit.searchWindow = 30;
cfg.orbit.windowGrowth = 1;
cfg.orbit.minCrossings = 4;
cfg.orbit.transientFraction = 0.2;
cfg.orbit.samplesPerCycle = 200;
cfg.orbit.extractNumPoints = 201;
model = struct('rhs', @circadian_linear_rhs);

orbit = circadian.find_orbit(model, [], cfg, cfg.initialState);

verifyTrue(testCase, orbit.success, 'Expected circadian find_orbit to return a periodic orbit.');
verifyEqual(testCase, numel(orbit.obs), numel(orbit.t));

features = evaluate_orbit_features(struct('t', orbit.t, 'y', orbit.y, 'period', orbit.period), ...
    @(t, y) y(:, 2) + y(:, 3), ...
    ["period", "observable_max", "observable_min", "observable_amplitude"], ...
    struct('extremaRefinement', false));
verifyEqual(testCase, orbit.obs, features.observable.series(:, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, orbit.amplitude, features.observable.amplitude(1), 'AbsTol', 2e-10);
verifyEqual(testCase, orbit.yMax, features.observable.max(1), 'AbsTol', 2e-10);
verifyEqual(testCase, orbit.yMin, features.observable.min(1), 'AbsTol', 2e-10);
end

function cfg = make_flexmod_cfg()
cfg = struct();
cfg.initialState = [1; 0];
cfg.orbit = struct( ...
    'searchWindow', 30, ...
    'windowGrowth', 1, ...
    'minCrossings', 4, ...
    'transientFraction', 0.2, ...
    'stateInitPhaseFraction', 0.25, ...
    'refinePeriods', 4, ...
    'refineTransientPeriods', 1, ...
    'samplesPerCycle', 200, ...
    'extractNumPoints', 201, ...
    'solverTol', struct('RelTol', 1e-6, 'AbsTol', 1e-9));
end

function dydt = harmonic_rhs(y, ~)
y = y(:);
dydt = [y(2); -y(1)];
end

function dydt = circadian_linear_rhs(y, ~)
y = y(:);
dydt = [y(2); -y(1); -y(1)];
end
