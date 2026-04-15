function tests = test_exported_views_for_postprocessing
%TEST_EXPORTED_VIEWS_FOR_POSTPROCESSING Exported solver/derived views support downstream callers.

    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(testDir);
    addpath(rootDir, '-begin');
    addpath(fullfile(rootDir, 'flexible_modulators'), '-begin');
    addpath(fullfile(rootDir, 'Circadian'), '-begin');
    testCase.TestData.rootDir = rootDir;
end

function testFlexmodOrbitReconstructionAcceptsDerivedView(testCase)
    task = make_reference_task();
    solverView = task.exportSolverView();
    derived = task.exportDerivedView();

    orbit = flexmod.orbit_from_state(derived, solverView.params);

    verifyEqual(testCase, orbit.period, derived.period, 'AbsTol', 1e-12);
    verifyEqual(testCase, orbit.params, solverView.params, 'AbsTol', 1e-12);
    verifyEqual(testCase, orbit.t(1:end-1), derived.t, 'AbsTol', 1e-12);
    verifyEqual(testCase, orbit.y(1:end-1,:), derived.TS_var, 'AbsTol', 1e-12);
    verifyEqual(testCase, orbit.t(end), derived.period, 'AbsTol', 1e-12);
    verifyEqual(testCase, orbit.y(end,:), derived.TS_var(1,:), 'AbsTol', 1e-12);
end

function testCircadianInitialGuessAcceptsDerivedView(testCase)
    task = make_reference_task();
    derived = task.exportDerivedView();
    cfg = struct('orbit', struct('stateInitPhaseFraction', 0.25));

    y0 = circadian.state_initial_guess(derived, cfg);

    sampleIndex = 1 + floor(cfg.orbit.stateInitPhaseFraction * max(size(derived.TS_var, 1) - 1, 0));
    sampleIndex = min(max(sampleIndex, 1), size(derived.TS_var, 1));
    verifyEqual(testCase, y0, derived.TS_var(sampleIndex, :).', 'AbsTol', 1e-12);
end

function testCircadianBuildTaskInputsAcceptsSolverAndDerivedViews(testCase)
    cfg = circadian.default_config();
    model = circadian.build_model(cfg);
    period0 = 24;
    solverView = struct( ...
        'params', reshape(model.defaultParams, 1, []), ...
        'p_Psi', [period0 / (2 * pi); zeros(3, 1)]);
    derivedView = struct( ...
        'obsAmp', 0.6, ...
        'obsMax', 1.4);

    taskSpec = struct();
    taskSpec.goalOrder = {'obsAmp', 'period'};
    taskSpec.controlledParams = {model.paramNames{1}, model.paramNames{2}};
    taskSpec.goals = struct('obsAmp', 0.9, 'period', 28);
    taskSpec.maxStepScale = 1;

    [itemsPerturb, itemsControlled, maxStep, accuracy] = ...
        circadian.build_task_inputs(taskSpec, model, solverView, derivedView, cfg);

    verifyEqual(testCase, itemsPerturb(1).prop, 'obsAmp');
    verifyEqual(testCase, itemsPerturb(1).idx, 1);
    verifyEqual(testCase, itemsPerturb(1).target, taskSpec.goals.obsAmp, 'AbsTol', 1e-12);
    verifyEqual(testCase, itemsPerturb(2).prop, 'p_Psi');
    verifyEqual(testCase, itemsPerturb(2).idx, 1);
    verifyEqual(testCase, itemsPerturb(2).target, taskSpec.goals.period / (2 * pi), 'AbsTol', 1e-12);
    verifyEqual(testCase, itemsControlled, ...
        [circadian.parameter_index(model, model.paramNames{1}), circadian.parameter_index(model, model.paramNames{2})]);
    verifyGreaterThan(testCase, maxStep(1), 0);
    verifyGreaterThan(testCase, maxStep(2), 0);
    verifyEqual(testCase, accuracy(1), cfg.fmam.targetAccuracy.obsAmp, 'AbsTol', 1e-12);
    verifyEqual(testCase, accuracy(2), cfg.fmam.targetAccuracy.period, 'AbsTol', 1e-12);
end

function testConstructorAcceptsSolverViewWithEquivalentBehavior(testCase)
    task0 = make_reference_task();
    solverView0 = task0.exportSolverView();
    derived0 = task0.exportDerivedView();

    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    item = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView0.p_Psi(1));

    task1 = FMAM_ODE(sys, obs, solverView0, item, 1, 0.1, 1e-6, 'derivatives', derivatives);

    verifyEqual(testCase, task1.exportSolverView().params, solverView0.params, 'AbsTol', 1e-12);
    verifyEqual(testCase, task1.exportSolverView().p_Psi, solverView0.p_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, task1.exportSolverView().q_Psi, solverView0.q_Psi, 'AbsTol', 1e-12);
    verifyEqual(testCase, task1.exportDerivedView().TS_var, derived0.TS_var, 'AbsTol', 1e-12);
    verifyEqual(testCase, task1.items_per_curr, task0.items_per_curr, 'AbsTol', 1e-12);

    fit0 = task0.fit();
    fit1 = task1.fit();
    verifyEqual(testCase, fit1.converged, fit0.converged);
    verifyEqual(testCase, task1.exportSolverView().p_Psi, task0.exportSolverView().p_Psi, 'AbsTol', 1e-12);

    task0.step();
    task1.step();
    verifyEqual(testCase, task1.continuationStatus.completed, task0.continuationStatus.completed);
    verifyEqual(testCase, task1.continuationStatus.lambda, task0.continuationStatus.lambda, 'AbsTol', 1e-12);
    verifyEqual(testCase, task1.exportDerivedView().period, task0.exportDerivedView().period, 'AbsTol', 1e-12);
end

function testExportedViewsSupportLegacyReprocessingInputs(testCase)
    task = make_reference_task();
    solverView = task.exportSolverView();
    derived = task.exportDerivedView();

    y0 = derived.TS_var(1, :).';
    searchWindow = max(10 * derived.period, derived.t(end));

    verifyEqual(testCase, reshape(solverView.params, 1, []), 1, 'AbsTol', 1e-12);
    verifyEqual(testCase, y0, derived.TS_var(1, :).', 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, searchWindow, derived.period);
end

function testConstructorAssemblySampleCountSetsDiscretization(testCase)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);
    item = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1));

    task = FMAM_ODE(sys, obs, solverView, item, 1, 0.1, 1e-6, ...
        'derivatives', derivatives, 'assemblySampleCount', 128);

    verifyEqual(testCase, task.assemblySampleCount, 128);
    verifyEqual(testCase, task.discretization.assemblySampleCount, 128);
end

function testConstructorRejectsLegacyLconstOption(testCase)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);
    item = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1));

    verifyError(testCase, @() FMAM_ODE(sys, obs, solverView, item, 1, 0.1, 1e-6, ...
        'derivatives', derivatives, 'Lconst', 128), ...
        'FMAM_ODE:InvalidConstructorOption');
end

function testReconstructionSettingsDriveExportedAndRebuiltViews(testCase)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);
    item = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1));
    reconstruction = struct('timeResampleCount', 4096, 'phaseSampleCount', 64);

    task = FMAM_ODE(sys, obs, solverView, item, 1, 0.1, 1e-6, ...
        'derivatives', derivatives, 'reconstruction', reconstruction);
    derived = task.exportDerivedView();
    rebuilt = task.rebuildState();

    verifyEqual(testCase, size(derived.TS_var, 1), reconstruction.phaseSampleCount);
    verifyEqual(testCase, task.reconstruction.phaseSampleCount, reconstruction.phaseSampleCount);
    verifyEqual(testCase, rebuilt.discretization.reconstruction.phaseSampleCount, ...
        reconstruction.phaseSampleCount);
    verifyEqual(testCase, size(rebuilt.TS_var, 1), reconstruction.phaseSampleCount);
end

function task = make_reference_task()
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);
    item = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1));
    task = FMAM_ODE(sys, obs, solverView, item, 1, 0.1, 1e-6, 'derivatives', derivatives);
end

function solverView = make_reference_solver_view(obs, params, M, PV)
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    [obs, params, t, x] = canonicalize_trajectory_fixture(obs, params, t, x, M, PV);
    discretization = fmam_state_defaults.defaultDiscretization();
    extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs, params, t, x, M, PV, ...
        discretization, extremaSearch);
end

function [obs, params, t, x] = canonicalize_trajectory_fixture(obs, params, t, x, M, PV)
    if isempty(obs)
        obs = {};
    elseif iscell(obs)
        obs = reshape(obs, 1, []);
    end
    params = params(:).';
    t = t(:);
    [t, x] = fmam_state_ops.normalizePeriodicInputs(t, x);
    fmam_state_ops.validateTrajectoryInputs(obs, params, t, x, M, PV);
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
