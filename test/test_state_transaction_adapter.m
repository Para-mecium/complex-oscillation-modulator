function tests = test_state_transaction_adapter
%TEST_STATE_TRANSACTION_ADAPTER Contract tests for solverView export/load and explicit state rebuild APIs.

    tests = functiontests(localfunctions);
end

function testExportLoadRoundTripRestoresSolverView(testCase)
    task = make_reference_task();
    snapshot0 = task.exportSolverView();

    mutated = snapshot0;
    mutated.p_Psi = mutated.p_Psi + linspace(-0.12, -0.01, numel(mutated.p_Psi)).';
    mutated.q_Psi = mutated.q_Psi + linspace(0.01, -0.005, numel(mutated.q_Psi)).';

    task.loadSolverView(mutated);
    mutatedView = task.exportSolverView();
    verifyNotEqual(testCase, mutatedView.p_Psi, snapshot0.p_Psi);
    verifyNotEqual(testCase, mutatedView.q_Psi, snapshot0.q_Psi);

    task.loadSolverView(snapshot0);
    verify_solver_view(testCase, task.exportSolverView(), snapshot0);
end

function testLoadPreservesMatrixBlockShapes(testCase)
    task = make_reference_task();
    snapshot0 = task.exportSolverView();

    mutated = snapshot0;
    delta = reshape((1:numel(mutated.p_var)).' / 10, size(mutated.p_var));
    mutated.p_var = mutated.p_var + delta;

    task.loadSolverView(mutated);
    mutatedView = task.exportSolverView();
    verifyEqual(testCase, mutatedView.p_var, snapshot0.p_var + delta, 'AbsTol', 1e-12);

    task.loadSolverView(snapshot0);
    verify_solver_view(testCase, task.exportSolverView(), snapshot0);
end

function testRebuildStateMatchesExportedResults(testCase)
    task = make_reference_task();
    stat0 = task.rebuildState();
    stat1 = task.rebuildState();
    derived = task.exportDerivedView();
    view = task.exportSolverView();

    verifyFalse(testCase, stat0 == stat1);
    verifyEqual(testCase, stat0.params, view.params, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat0.p_Psi, view.p_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat0.varPhiMax, view.varPhiMax, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat0.TS_var, derived.TS_var, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat0.period, derived.period, 'AbsTol', 1e-12);
    verifyEqual(testCase, stat0.t, derived.t, 'AbsTol', 1e-12);
end

function testStateFromViewsMatchesCompatibilityAdapter(testCase)
    task = make_reference_task();
    view = task.exportSolverView();
    derived = task.exportDerivedView();

    rebuilt = state.fromViews(task.obs, view, derived, task.discretization, task.extremaSearch);
    compat = task.rebuildState();

    verifyEqual(testCase, rebuilt.params, compat.params, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.p_Psi, compat.p_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.q_Psi, compat.q_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.TS_var, compat.TS_var, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.TS_obs, compat.TS_obs, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.period, compat.period, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.varAmp, compat.varAmp, 'AbsTol', 1e-12);
end

function testStateFromSolverViewMatchesCompatibilityAdapter(testCase)
    task = make_reference_task();
    view = task.exportSolverView();

    rebuilt = state.fromSolverView(task.obs, view, task.discretization, task.extremaSearch);
    compat = task.rebuildState();

    verifyEqual(testCase, rebuilt.params, compat.params, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.p_Psi, compat.p_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.TS_var, compat.TS_var, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.period, compat.period, 'AbsTol', 1e-12);
end

function testStateFromSolverSnapshotMatchesCompatibilityAdapter(testCase)
    task = make_reference_task();
    compat = task.rebuildState();
    snapshot = fmam_state_ops.buildStateSnapshotFromViews( ...
        task.exportSolverView(), task.exportDerivedView());

    rebuilt = state.fromSolverSnapshot(task.obs, snapshot, task.discretization, task.extremaSearch);

    verifyEqual(testCase, rebuilt.params, compat.params, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.p_Psi, compat.p_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.TS_var, compat.TS_var, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.TS_obs, compat.TS_obs, 'AbsTol', 1e-12);
    verifyEqual(testCase, rebuilt.period, compat.period, 'AbsTol', 1e-12);
end

function testConstructorRejectsLegacyStateInput(testCase)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    legacyState = state(obs, 1, t, x, 3, PV);
    items_per = struct('prop', 'p_Psi', 'idx', 1, 'target', legacyState.p_Psi(1));

    verifyError(testCase, ...
        @() FMAM_ODE(sys, obs, legacyState, items_per, 1, 0.1, 1e-6, 'derivatives', derivatives), ...
        'FMAM_ODE:InvalidInitialSolverView');
end

function task = make_reference_task(initialInput)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    if nargin < 1 || isempty(initialInput)
        PV = struct('name', 'var', 'idx', 1);
        t = linspace(0, 2 * pi, 1001).';
        x = [cos(t), sin(t)];
        initialInput = fmam_state_ops.solverViewFromState(state(obs, 1, t, x, 3, PV));
    end

    items_per = struct('prop', 'p_Psi', 'idx', 1, 'target', initialInput.p_Psi(1));
    task = FMAM_ODE(sys, obs, initialInput, items_per, 1, 0.1, 1e-6, 'derivatives', derivatives);
end

function verify_solver_view(testCase, view, snapshot)
    fields = {'params', 'p_Psi', 'q_Psi', 'p_var', 'q_var', ...
        'varPhiMax', 'varPhiMin', 'obsPhiMax', 'obsPhiMin'};

    for i = 1:numel(fields)
        fieldName = fields{i};
        verifyEqual(testCase, view.(fieldName), snapshot.(fieldName), 'AbsTol', 1e-12);
    end
end

function sys = make_harmonic_system()
    sys = cell(1, 2);
    sys{1} = @rhs_x1;
    sys{2} = @rhs_x2;
end

function dx = rhs_x1(y, parameters)
    dx = -parameters(1) .* y(:, 2);
end

function dx = rhs_x2(y, parameters)
    dx = parameters(1) .* y(:, 1);
end
