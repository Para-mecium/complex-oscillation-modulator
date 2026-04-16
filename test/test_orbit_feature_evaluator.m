function tests = test_orbit_feature_evaluator
%TEST_ORBIT_FEATURE_EVALUATOR Regression checks for orbit feature extraction.

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

function testRejectsUnsupportedOrbitFields(testCase)
orbit = make_closed_orbit();
orbit.extra = 1;

verifyError(testCase, @() evaluate_orbit_features(orbit, [], [], struct()), ...
    'orbitfeature:UnsupportedOrbitField');
end

function testRejectsMissingOrbitFields(testCase)
orbit = rmfield(make_closed_orbit(), 'period');

verifyError(testCase, @() evaluate_orbit_features(orbit, [], [], struct()), ...
    'orbitfeature:MissingOrbitField');
end

function testRejectsMismatchedPeriod(testCase)
orbit = make_closed_orbit();
orbit.period = orbit.period + 0.1;

verifyError(testCase, @() evaluate_orbit_features(orbit, [], [], struct()), ...
    'orbitfeature:InvalidOrbitInput');
end

function testStateFeaturesNormalizeRepeatedEndpoint(testCase)
orbit = make_closed_orbit();

features = evaluate_orbit_features(orbit, [], [], struct());

verifyEqual(testCase, features.period, orbit.period, 'AbsTol', 1e-12);
verifySize(testCase, features.state.series, [size(orbit.y, 1) - 1, size(orbit.y, 2)]);
verifyEqual(testCase, features.state.max, [1, 1], 'AbsTol', 2e-4);
verifyEqual(testCase, features.state.min, [-1, -1], 'AbsTol', 2e-4);
verifyEqual(testCase, features.state.amplitude, [1, 1], 'AbsTol', 2e-4);
end

function testRejectsComplexOrbitValues(testCase)
orbit = make_closed_orbit();
orbit.y(1, 1) = 1 + 1i;

verifyError(testCase, @() evaluate_orbit_features(orbit, [], [], struct()), ...
    'orbitfeature:InvalidOrbitInput');
end

function testObservableFeaturesSupportFunctionHandles(testCase)
orbit = make_closed_orbit();
obsSpec = { ...
    @(t, y) y(:, 2) + 0.5 * y(:, 1), ...
    @(t, y) [y(:, 1) - y(:, 2), y(:, 1) + y(:, 2)]};
featureSpec = ["period", "observable_max", "observable_min", "observable_amplitude"];

features = evaluate_orbit_features(orbit, obsSpec, featureSpec, struct());

verifySize(testCase, features.observable.series, [size(orbit.y, 1) - 1, 3]);
verifyEqual(testCase, features.observable.max(1), sqrt(1.25), 'AbsTol', 2e-4);
verifyEqual(testCase, features.observable.min(1), -sqrt(1.25), 'AbsTol', 2e-4);
verifyEqual(testCase, features.observable.amplitude(1), sqrt(1.25), 'AbsTol', 2e-4);
verifyEqual(testCase, features.observable.amplitude(2:3), [sqrt(2), sqrt(2)], 'AbsTol', 2e-4);
end

function testObservableFeaturesRequireObservableSpec(testCase)
orbit = make_closed_orbit();

verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], "observable_amplitude", struct()), ...
    'orbitfeature:ObservableRequired');
end

function testRejectsComplexObservableOutput(testCase)
orbit = make_closed_orbit();

verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, @(t, y) y(:, 1) + 1i * y(:, 2), ...
        "observable_amplitude", struct()), ...
    'orbitfeature:InvalidObservableOutput');
end

function testRejectsInvalidRefinementFactor(testCase)
orbit = make_closed_orbit();

verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], [], struct('refinementFactor', 1.5)), ...
    'orbitfeature:InvalidOptions');
verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], [], struct('refinementFactor', 1)), ...
    'orbitfeature:InvalidOptions');
end

function testRejectsInvalidRefinementPointCount(testCase)
orbit = make_closed_orbit();

verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], [], struct('refinementPointCount', 2)), ...
    'orbitfeature:InvalidOptions');
verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], [], struct('refinementPointCount', 4)), ...
    'orbitfeature:InvalidOptions');
verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], [], struct('refinementPointCount', 1.5)), ...
    'orbitfeature:InvalidOptions');
end

function testRejectsInvalidRefinementMethod(testCase)
orbit = make_closed_orbit();

verifyError(testCase, ...
    @() evaluate_orbit_features(orbit, [], [], struct('refinementMethod', "linear")), ...
    'orbitfeature:InvalidOptions');
end

function testSplineRefinementMethodIsSupported(testCase)
orbit = make_closed_orbit();

features = evaluate_orbit_features(orbit, [], ...
    ["state_max", "state_min", "state_amplitude"], ...
    struct('refinementMethod', "spline", 'refinementPointCount', 5));

verifyTrue(testCase, all(features.diagnostics.state.maxRefined));
verifyTrue(testCase, all(features.diagnostics.state.minRefined));
verifyEqual(testCase, features.state.amplitude, [1, 1], 'AbsTol', 5e-3);
end

function testQuadraticRefinementMethodIsSupported(testCase)
orbit = make_closed_orbit();

features = evaluate_orbit_features(orbit, [], ...
    ["state_max", "state_min", "state_amplitude"], ...
    struct('refinementMethod', "quadratic", 'refinementPointCount', 5));

verifyTrue(testCase, all(features.diagnostics.state.maxRefined));
verifyTrue(testCase, all(features.diagnostics.state.minRefined));
verifyEqual(testCase, features.state.amplitude, [1, 1], 'AbsTol', 5e-3);
end

function testStateFeatureAccuracyAgainstAnalyticSinusoids(testCase)
period = 2 * pi;
t = linspace(0, period, 65).';

a1 = 1.3;
b1 = 0.4;
c1 = 0.2;
a2 = 0.6;
b2 = -1.1;
c2 = -0.7;

orbit = struct( ...
    't', t, ...
    'y', [c1 + a1 * sin(t) + b1 * cos(t), ...
          c2 + a2 * sin(t) + b2 * cos(t)], ...
    'period', period);

features = evaluate_orbit_features(orbit, [], ...
    ["period", "state_max", "state_min", "state_amplitude"], struct());

ampTrue = [hypot(a1, b1), hypot(a2, b2)];
maxTrue = [c1, c2] + ampTrue;
minTrue = [c1, c2] - ampTrue;

verifyEqual(testCase, features.period, period, 'AbsTol', 1e-12);
verifyEqual(testCase, features.state.max, maxTrue, 'AbsTol', 5e-4);
verifyEqual(testCase, features.state.min, minTrue, 'AbsTol', 5e-4);
verifyEqual(testCase, features.state.amplitude, ampTrue, 'AbsTol', 5e-4);
end

function testExtremaRefinementHasBoundedShiftedPeakError(testCase)
period = 2 * pi;
t = linspace(0, period, 17).';
phaseShift = 0.37;
offset = -0.4;
amplitudeTrue = 1.7;
maxTrue = offset + amplitudeTrue;
minTrue = offset - amplitudeTrue;
orbit = struct( ...
    't', t, ...
    'y', offset + amplitudeTrue * sin(t + phaseShift), ...
    'period', period);

refined = evaluate_orbit_features(orbit, [], ...
    ["state_max", "state_min", "state_amplitude"], struct());

refinedAmpError = abs(refined.state.amplitude - amplitudeTrue);
refinedMaxError = abs(refined.state.max - maxTrue);
refinedMinError = abs(refined.state.min - minTrue);

verifyLessThan(testCase, refinedAmpError, 5e-3);
verifyLessThan(testCase, refinedMaxError, 5e-3);
verifyLessThan(testCase, refinedMinError, 5e-3);
verifyTrue(testCase, refined.diagnostics.state.maxRefined);
verifyTrue(testCase, refined.diagnostics.state.minRefined);
end

function testRefinementFallbackIsReported(testCase)
orbit = struct( ...
    't', [0; 1; 2], ...
    'y', [0; 1; 0], ...
    'period', 2);

features = evaluate_orbit_features(orbit, [], [], struct());

verifyTrue(testCase, all(features.diagnostics.state.maxFallback));
verifyTrue(testCase, all(features.diagnostics.state.minFallback));
verifyEqual(testCase, features.state.amplitude, 0.5, 'AbsTol', 1e-12);
end

function testRefinementFallsBackWhenPointCountExceedsSamples(testCase)
orbit = struct( ...
    't', [0; 1; 2; 3], ...
    'y', [0; 1; -1; 0], ...
    'period', 3);

features = evaluate_orbit_features(orbit, [], [], struct('refinementPointCount', 5));

verifyTrue(testCase, all(features.diagnostics.state.maxFallback));
verifyTrue(testCase, all(features.diagnostics.state.minFallback));
verifyEqual(testCase, features.diagnostics.state.maxReason, "insufficient_local_points");
verifyEqual(testCase, features.diagnostics.state.minReason, "insufficient_local_points");
end

function orbit = make_closed_orbit()
period = 2 * pi;
t = linspace(0, period, 257).';
orbit = struct( ...
    't', t, ...
    'y', [sin(t), cos(t)], ...
    'period', period);
end
