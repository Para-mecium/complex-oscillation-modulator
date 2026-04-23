function tests = test_state_regressions
%TEST_STATE_REGRESSIONS Regression coverage for state preprocessing.

    tests = functiontests(localfunctions);
end

function testTimeOriginShiftDoesNotChangeRecoveredOrbit(testCase)
    T = 10;
    M = 5;
    obs = {};
    params = 1;
    PV = struct('name', 'var', 'idx', 1);

    tBase = linspace(0, T, 1002).';
    tBase(end) = [];
    x = [cos(2 * pi * tBase / T), sin(2 * pi * tBase / T)];

    state0 = state(obs, params, tBase, x, M, PV);
    stateShifted = state(obs, params, tBase + 5, x, M, PV);

    verifyEqual(testCase, stateShifted.period, state0.period, 'AbsTol', 5e-2);
    verifyEqual(testCase, stateShifted.p_var, state0.p_var, 'AbsTol', 5e-3);
    verifyEqual(testCase, stateShifted.q_var, state0.q_var, 'AbsTol', 5e-3);
end

function testClosedEndpointSamplesAreAccepted(testCase)
    T = 2 * pi;
    obs = {};
    params = 1;
    PV = struct('name', 'var', 'idx', 1);
    t = linspace(0, T, 1001).';
    x = [cos(t), sin(t)];

    stat = state(obs, params, t, x, 5, PV);

    verifyEqual(testCase, stat.period, T, 'AbsTol', 5e-2);
    verifyEqual(testCase, stat.p_var(1, 1), 0, 'AbsTol', 5e-2);
    verifyEqual(testCase, stat.p_var(2, 1), 1, 'AbsTol', 5e-2);
    verifyEqual(testCase, stat.q_var(:, 1), zeros(size(stat.q_var(:, 1))), 'AbsTol', 5e-2);
end

function testLinearTimePhaseModeKeepsConstantPsiAndNativeFourierContent(testCase)
    T = 10;
    M = 5;
    obs = {};
    params = 1;
    PV = struct('name', 'var', 'idx', 1);
    discretization = fmam_state_defaults.defaultDiscretization();
    discretization.reconstruction.phaseMode = 'linearTime';
    discretization.reconstruction.phaseSampleCount = 256;
    discretization.reconstruction.timeResampleCount = 512;

    t = linspace(0, T, 1001).';
    theta = 2 * pi * t / T;
    x = [cos(theta) + 0.2 * cos(2 * theta), sin(theta)];

    stat = state(obs, params, t, x, M, PV, discretization);

    verifyEqual(testCase, stat.period, T, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat.p_Psi(1), T / (2 * pi), 'AbsTol', 1e-12);
    verifyEqual(testCase, stat.p_Psi(2:end), zeros(M - 1, 1), 'AbsTol', 1e-12);
    verifyEqual(testCase, stat.q_Psi, zeros(M - 1, 1), 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, abs(stat.p_var(3, 1)), 0.15);
    verifyLessThan(testCase, abs(stat.p_var(2, 1) - 1), 1e-10);
end

function testLinearTimePhaseModeCanBePassedByNameValue(testCase)
    T = 2 * pi;
    M = 4;
    obs = {};
    params = 1;
    PV = struct('name', 'var', 'idx', 1);
    discretization = fmam_state_defaults.defaultDiscretization();
    discretization.reconstruction.phaseMode = 'linearTime';
    t = linspace(0, T, 101).';
    x = [cos(t), sin(t)];

    stat = state(obs, params, t, x, M, PV, ...
        'discretization', discretization);

    verifyEqual(testCase, stat.p_Psi(1), 1, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat.q_Psi, zeros(M - 1, 1), 'AbsTol', 1e-12);
end

function testNonMonotoneTimeMapEmitsNamedWarning(testCase)
    stat = make_reference_state();
    snapshot = make_solver_snapshot(stat);
    snapshot.p_Psi = zeros(size(snapshot.p_Psi));
    snapshot.q_Psi = zeros(size(snapshot.q_Psi));
    stat = state.fromSolverSnapshot(stat.obs, snapshot, stat.discretization, stat.extremaSearch);

    verifyWarning(testCase, @() stat.assertTimeMapInvariant(), 'state:NonMonotoneTimeMap');

    statForTimeGrid = state.fromSolverSnapshot( ...
        stat.obs, snapshot, stat.discretization, stat.extremaSearch);
    verifyWarning(testCase, @() statForTimeGrid.t, 'state:NonMonotoneTimeMap');
end

function testNonMonotoneTimeMapCheckCanBeDisabled(testCase)
    stat = make_reference_state();
    snapshot = make_solver_snapshot(stat);
    snapshot.p_Psi = zeros(size(snapshot.p_Psi));
    snapshot.q_Psi = zeros(size(snapshot.q_Psi));
    stat = state.fromSolverSnapshot(stat.obs, snapshot, stat.discretization, stat.extremaSearch);
    stat.checkPsiNonnegative = false;

    verifyWarningFree(testCase, @() stat.assertTimeMapInvariant());
    verifyWarningFree(testCase, @() stat.t);
end

function testTimeGridUsesExactIntervalIntegrals(testCase)
    stat = make_reference_state();
    phaseSampleCount = stat.discretization.reconstruction.phaseSampleCount;
    dt = fmam_state_ops.timeIncrementsFromCoefficients(stat.p_Psi, stat.q_Psi, phaseSampleCount);
    tGrid = stat.t;
    tClosed = [tGrid; stat.period];

    verifyGreaterThan(testCase, min(dt), 0);
    verifyEqual(testCase, diff(tClosed), dt, 'AbsTol', 1e-10);
end

function testDiscretizationSetterMergesPartialAndWholeStructUpdates(testCase)
    stat = make_reference_state();
    original = stat.discretization;

    stat.discretization = struct('reconstruction', struct('phaseSampleCount', 64));
    verifyEqual(testCase, stat.discretization.reconstruction.phaseSampleCount, 64);
    verifyEqual(testCase, stat.discretization.reconstruction.timeResampleCount, ...
        original.reconstruction.timeResampleCount);
    verifyEqual(testCase, numel(stat.phi), 64);
    verifyEqual(testCase, size(stat.TS_var, 1), 64);

    stat.discretization = struct( ...
        'assemblySampleCount', 256, ...
        'reconstruction', struct('timeResampleCount', 4096, 'phaseSampleCount', 96));
    verifyEqual(testCase, stat.discretization.assemblySampleCount, 256);
    verifyEqual(testCase, stat.discretization.reconstruction.timeResampleCount, 4096);
    verifyEqual(testCase, stat.discretization.reconstruction.phaseSampleCount, 96);
    verifyEqual(testCase, numel(stat.phi), 96);
end

function testExtremaSearchSetterMergesPartialAndWholeStructUpdates(testCase)
    stat = make_reference_state();
    original = stat.extremaSearch;

    stat.extremaSearch.maxRefinementRounds = 7;
    verifyEqual(testCase, stat.extremaSearch.maxRefinementRounds, 7);
    verifyEqual(testCase, stat.extremaSearch.extremaResidualTolerance, ...
        original.extremaResidualTolerance);

    stat.extremaSearch = struct( ...
        'maxRefinementRounds', 4, ...
        'extremaResidualTolerance', 1e-5);
    verifyEqual(testCase, stat.extremaSearch.maxRefinementRounds, 4);
    verifyEqual(testCase, stat.extremaSearch.extremaResidualTolerance, 1e-5);
end

function testValidateTrajectoryInputsRejectsNonMonotoneTimeGrid(testCase)
    t = [0; 0.5; 0.25];
    x = [1 0; 0 1; -1 0];

    didFail = false;
    try
        fmam_state_ops.validateTrajectoryInputs({}, 1, t, x, 3, struct('name', 'var', 'idx', 1));
    catch ME
        didFail = true;
        verifySubstring(testCase, ME.message, 'strictly increasing');
    end

    verifyTrue(testCase, didFail);
end

function stat = make_reference_state()
    t = linspace(0, 2 * pi, 1001).';
    t(end) = [];
    x = [cos(t), sin(t)];
    stat = state({}, 1, t, x, 5, struct('name', 'var', 'idx', 1));
end

function snapshot = make_solver_snapshot(stat)
    snapshot = fmam_state_ops.buildStateSnapshotFromViews( ...
        fmam_state_ops.solverViewFromState(stat), ...
        fmam_state_ops.derivedViewFromState(stat));
end
