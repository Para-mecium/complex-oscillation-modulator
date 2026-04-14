function tests = test_data_driven_regression
%TEST_DATA_DRIVEN_REGRESSION Regression coverage for the data-driven FMAM setup.

    tests = functiontests(localfunctions);
end

function testInitialDataDrivenFitIsNotBlockedByIntermediatePsiChecks(testCase)
    [task, initialSolverView] = make_data_driven_task();

    result = task.fit();

    verifyTrue(testCase, result.converged);
    verifyLessThanOrEqual(testCase, max(result.finalError), task.errBound);
    verifyGreaterThanOrEqual(testCase, result.iterations, 1);
    verifyEqual(testCase, size(initialSolverView.p_Psi), [task.truncationOrder, 1]);
    verifyEqual(testCase, size(initialSolverView.q_Psi), [task.truncationOrder - 1, 1]);
end

function testDataDrivenFitMatchesLegacyBehaviorWithinTolerance(testCase)
    baseline = legacy_data_driven_baseline();
    [task, ~] = make_data_driven_task();

    result = task.fit();
    solverView = task.exportSolverView();
    derived = task.exportDerivedView();
    targetProperties = current_target_properties(solverView, derived);
    controlledParameters = solverView.params([1 4 5 6]);

    verifyTrue(testCase, result.converged);
    for i = 1:numel(baseline.targetProperties)
        absTol = max(5e-6, 5e-3 * abs(baseline.targetProperties(i)));
        verifyEqual(testCase, targetProperties(i), baseline.targetProperties(i), 'AbsTol', absTol);
    end

    for i = 1:numel(baseline.controlledParameters)
        absTol = max(1e-4, 1e-2 * abs(baseline.controlledParameters(i)));
        verifyEqual(testCase, controlledParameters(i), baseline.controlledParameters(i), 'AbsTol', absTol);
    end

    verifyLessThanOrEqual(testCase, max(result.finalError), ...
        max(task.errBound, 10 * baseline.legacyFinalMaxError + 1e-12));
end

function testDataDrivenStepUsesInitialLambdaStepWhenNoMaxStepCap(testCase)
    [task, ~] = make_data_driven_task();
    task.needLog = true;

    output = evalc('task.step();');

    verifyNotEmpty(testCase, task.logs);
    verifyNotEmpty(testCase, strfind(output, 'initial dlambda=1.000000e-01')); %#ok<STRIFCND>
    verifyNotEmpty(testCase, strfind(output, 'lambda 0.000000 -> 0.100000')); %#ok<STRIFCND>
end

function [task, initialSolverView] = make_data_driven_task()
    testDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(testDir);
    circuitDir = fullfile(rootDir, 'RLT_circuit');
    startDir = pwd;
    cleanupObj = onCleanup(@() cd(startDir)); %#ok<NASGU>

    cd(circuitDir);
    load("initData_circuit.mat");
    System;
    obs = [];
    load("initData_ODE.mat");

    params = Parameters;
    t = TS{1};
    TS_var = TS{2};
    M = 50;

    cd(rootDir);
    derivatives = build_symbolic_derivatives(sys, obs, numel(params));
    PV = struct('name', 'var', 'idx', 1);
    [obs, params, t, TS_var] = canonicalize_trajectory_fixture(obs, params, t, TS_var, M, PV);
    discretization = state.defaultDiscretization();
    extremaSearch = state.defaultExtremaSearch();
    initialSolverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs, params, t, TS_var, M, PV, ...
        discretization, extremaSearch);

    items_per = struct;
    items_per(1).prop = 'p_Psi';
    items_per(1).idx = 1;
    items_per(1).target = period / (2 * pi);

    items_per(2).prop = 'varAmp';
    items_per(2).idx = 1;
    items_per(2).target = varAmp(1);

    items_per(3).prop = 'varAmp';
    items_per(3).idx = 2;
    items_per(3).target = varAmp(2);

    items_per(4).prop = 'varAmp';
    items_per(4).idx = 3;
    items_per(4).target = varAmp(3);

    items_controlled = [1 4 5 6];
    errBound = 1e-6;
    continuationOptions = struct('initialLambdaStep', 0.1);

    task = FMAM_ODE(sys, obs, initialSolverView, items_per, items_controlled, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    task.isPsiUpdated = true;
end

function [obs, params, t, TS_var] = canonicalize_trajectory_fixture(obs, params, t, TS_var, M, PV)
    if isempty(obs)
        obs = {};
    elseif iscell(obs)
        obs = reshape(obs, 1, []);
    end
    params = params(:).';
    t = t(:);
    [t, TS_var] = fmam_state_ops.normalizePeriodicInputs(t, TS_var);
    fmam_state_ops.validateTrajectoryInputs(obs, params, t, TS_var, M, PV);
end

function values = current_target_properties(solverView, derived)
    values = zeros(1, 4);
    values(1) = solverView.p_Psi(1);
    values(2:4) = reshape(derived.varAmp(1:3), 1, []);
end

function baseline = legacy_data_driven_baseline()
    baseline = struct();
    baseline.targetProperties = [ ...
        0.712202534934838, ...
        0.6876636248430121, ...
        0.6876241174877468, ...
        0.6902592099613101];
    baseline.controlledParameters = [ ...
        500.1247867482443, ...
        0.9983394718571603, ...
        1.000173518651478, ...
        1.003891952414744];
    baseline.legacyFinalMaxError = 4.867217739956686e-13;
end
