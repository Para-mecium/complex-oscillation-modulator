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

function testNonMonotoneTimeMapEmitsNamedWarning(testCase)
    stat = make_reference_state();
    snapshot = stat.snapshotSolverState();
    snapshot.p_Psi = zeros(size(snapshot.p_Psi));
    snapshot.q_Psi = zeros(size(snapshot.q_Psi));
    stat.restoreSolverState(snapshot);

    verifyWarning(testCase, @() stat.assertTimeMapInvariant(), 'state:NonMonotoneTimeMap');

    statForTimeGrid = make_reference_state();
    statForTimeGrid.restoreSolverState(snapshot);
    verifyWarning(testCase, @() statForTimeGrid.t, 'state:NonMonotoneTimeMap');
end

function testNonMonotoneTimeMapCheckCanBeDisabled(testCase)
    stat = make_reference_state();
    snapshot = stat.snapshotSolverState();
    snapshot.p_Psi = zeros(size(snapshot.p_Psi));
    snapshot.q_Psi = zeros(size(snapshot.q_Psi));
    stat.restoreSolverState(snapshot);
    stat.checkPsiNonnegative = false;

    verifyWarningFree(testCase, @() stat.assertTimeMapInvariant());
    verifyWarningFree(testCase, @() stat.t);
end

function testTimeGridUsesExactIntervalIntegrals(testCase)
    stat = make_reference_state();
    dt = fmam_state_ops.timeIncrementsFromCoefficients(stat.p_Psi, stat.q_Psi, stat.LphiConst);
    tGrid = stat.t;
    tClosed = [tGrid; stat.period];

    verifyGreaterThan(testCase, min(dt), 0);
    verifyEqual(testCase, diff(tClosed), dt, 'AbsTol', 1e-10);
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
