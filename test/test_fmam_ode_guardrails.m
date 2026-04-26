function tests = test_fmam_ode_guardrails
%TEST_FMAM_ODE_GUARDRAILS Regression guardrails for FMAM_ODE construction and target semantics.

    tests = functiontests(localfunctions);
end

function testConstructorNormalizesConfigurationAndFreezesInitialTargets(testCase)
    task = make_guardrail_task();

    verifyEqual(testCase, size(task.items_perturb), [1 5]);
    verifyEqual(testCase, size(task.items_controlled), [1 5]);
    verifyEqual(testCase, task.items_controlled, [5 4 3 2 1]);
    verifyEqual(testCase, task.items_perturb(1).prop, 'params');
    verifyEqual(testCase, task.items_perturb(2).prop, 'p_Psi');
    verifyEqual(testCase, task.items_perturb(3).prop, 'varAmp');
    verifyEqual(testCase, task.items_perturb(4).prop, 'obsAmp');
    verifyEqual(testCase, task.items_perturb(5).prop, 'varPhase');
    verifyEqual(testCase, task.maxstepsize, 0.25 * ones(1, 5), 'AbsTol', 1e-12);
    verifyEqual(testCase, task.accuracy, 1e-6 * ones(1, 5), 'AbsTol', 1e-12);
    verifyEqual(testCase, task.continuationOptions.predictorMode, 'hermite');
    verifyEqual(testCase, task.continuationOptions.quadraticStepRatioBounds, [0.75 1.5]);
    verifyEqual(testCase, task.continuationOptions.initialLambdaStep, 0.2, 'AbsTol', 1e-12);

    expectedTargets = current_guardrail_target_values(task);
    verifyEqual(testCase, task.items_per_curr, expectedTargets, 'AbsTol', 1e-10);

    frozenTargets = task.items_per_curr;
    mutate_solver_view(task, @(view) with_field(view, 'params', replace_at(view.params, 2, view.params(2) + 3.5)));
    mutate_solver_view(task, @(view) with_field(view, 'p_Psi', replace_at(view.p_Psi, 1, view.p_Psi(1) + 0.25)));
    verifyEqual(testCase, task.items_per_curr, frozenTargets, 'AbsTol', 1e-12);
end

function testWrappedPeriodicTargetsUseShortestContinuationPath(testCase)
    task = make_wrapped_target_task();
    view = task.exportSolverView();
    initialTargets = [view.varPhiMax(1), view.obsPhiMax(1)];

    task.autostepsize();
    verifyEqual(testCase, task.items_per_curr, initialTargets, 'AbsTol', 1e-10);
    task.perturb();

    verifyEqual(testCase, task.items_per_curr, initialTargets + [-0.1 -0.1], ...
        'AbsTol', 1e-10);
end

function testConstructorRejectsUnknownNameValueOption(testCase)
    sys = make_harmonic_system();
    obs = {@observable_x1};
    params = ones(1, 5);
    derivatives = build_symbolic_derivatives(sys, obs, numel(params));
    [solverView, ~] = make_guardrail_views(obs, params);
    items_per = struct('prop', 'params', 'idx', 2, 'target', solverView.params(2));

    verifyError(testCase, ...
        @() FMAM_ODE(sys, obs, solverView, items_per, 1, 0.25, 1e-6, ...
            'derivatives', derivatives, ...
            'unknownOption', 1), ...
	        'FMAM_ODE:InvalidConstructorOption');
end

function testConstructorUsesPsiUpdateModeOptionOnly(testCase)
    task = make_guardrail_task('psiUpdateMode', true);

    verifyTrue(testCase, task.psiUpdateMode);
    verifyError(testCase, @() make_guardrail_task('isPsiUpdated', true), ...
        'FMAM_ODE:InvalidConstructorOption');
end

function testConstructorRejectsUnknownContinuationOptionField(testCase)
    verifyError(testCase, ...
        @() make_guardrail_task('continuationOptions', struct('predictorMod', 'constant')), ...
        'FMAM_ODE:UnknownContinuationOption');
end

function testConstructorRejectsUnknownNewtonOptionField(testCase)
    verifyError(testCase, ...
        @() make_guardrail_task('newtonOptions', struct('maxIteratoins', 1)), ...
        'FMAM_ODE:UnknownNewtonOption');
end

function testConstructorAcceptsNewtonStabilizationOptions(testCase)
    task = make_guardrail_task('newtonOptions', struct( ...
        'linearSystemScaling', 'none', ...
        'requireDescent', false, ...
        'acceptIncreaseTolerance', 0.25));

    verifyEqual(testCase, task.newtonOptions.linearSystemScaling, 'none');
    verifyFalse(testCase, task.newtonOptions.requireDescent);
    verifyEqual(testCase, task.newtonOptions.acceptIncreaseTolerance, 0.25, ...
        'AbsTol', 0);
end

function testConstructorRejectsInvalidNewtonScalingOption(testCase)
    verifyError(testCase, ...
        @() make_guardrail_task('newtonOptions', struct('linearSystemScaling', 'column')), ...
        'FMAM_ODE:InvalidNewtonOption');
end

function testRuntimeKnobsDoNotChangeFrozenTargetSnapshot(testCase)
    task = make_guardrail_task();
    frozenTargets = task.items_per_curr;

    task.needLog = ~task.needLog;
    task.verbose = ~task.verbose;
    task.errBound = task.errBound / 10;

    verifyEqual(testCase, task.items_per_curr, frozenTargets, 'AbsTol', 1e-12);

    mutate_solver_view(task, @(view) with_field(view, 'params', replace_at(view.params, 2, view.params(2) + 3.5)));
    mutate_solver_view(task, @(view) with_field(view, 'p_Psi', replace_at(view.p_Psi, 1, view.p_Psi(1) + 0.25)));
    verifyEqual(testCase, task.items_per_curr, frozenTargets, 'AbsTol', 1e-12);
end

function testAutostepsizeUsesFrozenObservableExtremaView(testCase)
    [task,armMutation] = make_mutating_obs_amp_task();
    frozenTarget = task.items_per_curr;

    armMutation(true);
    task.autostepsize();

    verifyEqual(testCase, task.items_per_curr, frozenTarget, 'AbsTol', 1e-12);
end

function testResidualUsesFrozenTargetViewWhenObservableMutatesState(testCase)
    [task,armMutation] = make_mutating_obs_amp_task();
    frozenTarget = task.items_per_curr;

    armMutation(true);
    residual = task.res();

    verifyEqual(testCase, residual(4), 0, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.items_per_curr, frozenTarget, 'AbsTol', 1e-12);
end

function testResidualIncludesObservablePrimaryGaugeViolation(testCase)
    task = make_observable_pv_guardrail_task();

    baselineResidual = task.res();
    baselineGaugeResidual = baselineResidual(end);

    mutate_solver_view(task, @(view) with_field(view, 'p_var', ...
        bump_matrix_entry(view.p_var, 3, 1, 5e-2)));
    residual = task.res();

    verifyGreaterThan(testCase, residual(end), baselineGaugeResidual + 1e-3);
end

function testPsiUpdateModeAssignmentsCaptureReferences(testCase)
    task = make_guardrail_task();

    verifyFalse(testCase, task.psiUpdateMode);

    task.psiUpdateMode = true;
    task.refreshPsiModeReferences();
    verifyTrue(testCase, task.psiUpdateMode);

    task.psiUpdateMode = false;
    task.refreshPsiModeReferences();
    verifyFalse(testCase, task.psiUpdateMode);
    verifyEqual(testCase, psi_reference_residual(task), ...
        zeros(max(0,2 * task.truncationOrder - 2), 1), 'AbsTol', 1e-12);
end

function testPsiNonnegativeCheckDefaultsOnAndCanBeConfigured(testCase)
    task = make_guardrail_task();

    verifyTrue(testCase, task.checkPsiNonnegative);

    task.checkPsiNonnegative = false;
    verifyFalse(testCase, task.checkPsiNonnegative);
    verifyFalse(testCase, task.rebuildState().checkPsiNonnegative);

    task.checkPsiNonnegative = true;
    verifyTrue(testCase, task.checkPsiNonnegative);

    configuredTask = make_guardrail_task('checkPsiNonnegative', false);
    verifyFalse(testCase, configuredTask.checkPsiNonnegative);
end

function testDisablePsiUpdateCapturesPsiReference(testCase)
    task = make_guardrail_task();
    delta_p = [0; 0.125; -0.05];
    delta_q = [0.05; -0.125];

    task.psiUpdateMode = true;
    task.refreshPsiModeReferences();
    mutate_solver_view(task, @(view) apply_psi_delta(view, delta_p, delta_q));

    task.psiUpdateMode = false;
    task.refreshPsiModeReferences();
    verifyFalse(testCase, task.psiUpdateMode);
    verifyEqual(testCase, psi_reference_residual(task), ...
        zeros(max(0,2 * task.truncationOrder - 2), 1), 'AbsTol', 1e-12);

    task.psiUpdateMode = false;
    task.refreshPsiModeReferences();
    verifyEqual(testCase, psi_reference_residual(task), ...
        zeros(max(0,2 * task.truncationOrder - 2), 1), 'AbsTol', 1e-12);
end

function testLoadSolverViewRefreshesFrozenReferencesByDefault(testCase)
    task = make_guardrail_task();
    delta_p = [0; 0.125; -0.05];
    delta_q = [0.05; -0.125];
    view = apply_psi_delta(task.exportSolverView(), delta_p, delta_q);

    task.loadSolverView(view);

    verifyEqual(testCase, psi_reference_residual(task), ...
        zeros(max(0,2 * task.truncationOrder - 2), 1), 'AbsTol', 1e-12);
end

function testLoadSolverViewCanPreserveFrozenReferences(testCase)
    task = make_guardrail_task();
    delta_p = [0; 0.125; -0.05];
    delta_q = [0.05; -0.125];
    view = apply_psi_delta(task.exportSolverView(), delta_p, delta_q);

    task.loadSolverView(view, struct('refreshPsiModeReferences', false));

    verifyEqual(testCase, psi_reference_residual(task), ...
        -[delta_p(2:end); delta_q(:)], 'AbsTol', 1e-12);
end

function testFrozenPhaseResidualParticipatesInConvergence(testCase)
    task = make_guardrail_task();
    task.errBound = inf;
    task.oneIter();

    mutate_solver_view(task, @(view) with_field(view, 'p_var', ...
        bump_matrix_entry(view.p_var, 3, 1, 5e-2)));
    task.errBound = 1e-12;
    result = task.oneIter();

    verifyFalse(testCase, result.converged);
    verifyGreaterThan(testCase, result.scalarError, task.errBound);
end

function testStructuralPropertiesCannotBeReboundAfterConstruction(testCase)
    task = make_guardrail_task();

    verifyPropertyWriteFails(testCase, @() assign_continuation_options(task));
    verifyPropertyWriteFails(testCase, @() assign_items_perturb(task));
end

function task = make_guardrail_task(varargin)
    sys = make_harmonic_system();
    obs = {@observable_x1};
    params = ones(1, 5);
    derivatives = build_symbolic_derivatives(sys, obs, numel(params));
    [solverView, derived] = make_guardrail_views(obs, params);

    items_per(1) = struct('prop', 'PARAMS', 'idx', 2, 'target', solverView.params(2));
    items_per(2) = struct('prop', 'p_psi', 'idx', 1, 'target', solverView.p_Psi(1));
    items_per(3) = struct('prop', 'VARAMP', 'idx', 1, 'target', derived.varAmp(1));
    items_per(4) = struct('prop', 'obsamp', 'idx', 1, 'target', derived.obsAmp(1));
    items_per(5) = struct('prop', 'VarPhase', 'idx', [1 2], 'target', derived.varPhase(1, 2));

    items_controlled = [5; 4; 3; 2; 1];
    continuationOptions = struct( ...
        'initialLambdaStep', 0.2, ...
        'predictorMode', "HERMITE", ...
        'quadraticStepRatioBounds', [0.75; 1.5]);

    task = FMAM_ODE(sys, obs, solverView, items_per, items_controlled, 0.25, 1e-6, ...
        'derivatives', derivatives, ...
        'continuationOptions', continuationOptions, varargin{:});
    task.needLog = false;
end

function task = make_wrapped_target_task()
    sys = make_harmonic_system();
    obs = {@observable_x1};
    params = ones(1, 5);
    derivatives = build_symbolic_derivatives(sys, obs, numel(params));
    [solverView, ~] = make_guardrail_views(obs, params);

    wrappedShift = 2 * pi - 0.1;
    items_per(1) = struct('prop', 'varPhiMax', 'idx', 1, ...
        'target', solverView.varPhiMax(1) + wrappedShift);
    items_per(2) = struct('prop', 'obsPhiMax', 'idx', 1, ...
        'target', solverView.obsPhiMax(1) + wrappedShift);

    task = FMAM_ODE(sys, obs, solverView, items_per, [1 2], 1, 1e-6, ...
        'derivatives', derivatives, ...
        'continuationOptions', struct('initialLambdaStep', 1, 'predictorMode', 'constant'));
    task.needLog = false;
end

function [task,armMutation] = make_mutating_obs_amp_task()
    sys = make_harmonic_system();
    baseObs = @observable_x1;
    obs = {[]};
    params = ones(1, 5);
    derivatives = build_symbolic_derivatives(sys, {baseObs}, numel(params));
    [solverView, derived] = make_guardrail_views({baseObs}, params);
    mutationProbe = make_guardrail_state({baseObs}, params);
    item = struct('prop', 'obsAmp', 'idx', 1, 'target', derived.obsAmp(1));
    armed = false;
    triggered = false;

    obs{1} = @mutating_observable_x1;
    task = FMAM_ODE(sys, obs, solverView, item, 1, 0.25, 1e-6, ...
        'derivatives', derivatives);
    task.needLog = false;

    armMutation = @set_mutation_arm;

    function z = mutating_observable_x1(y)
        if armed && ~triggered
            triggered = true;
            mutationProbe.add('obsPhiMin', 1, -mutationProbe.obsPhiMin(1));
        end
        z = y(:, 1);
    end

    function set_mutation_arm(flag)
        armed = logical(flag);
        if ~armed
            triggered = false;
        end
    end
end

function [solverView, derived] = make_guardrail_views(obs, params)
    PV = struct('name', 'var', 'idx', 1);
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    [obs, params, t, x] = canonicalize_trajectory_fixture(obs, params, t, x, 3, PV);
    discretization = fmam_state_defaults.defaultDiscretization();
    extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs, params, t, x, 3, PV, ...
        discretization, extremaSearch);
    derived = fmam_state_ops.buildDerivedView( ...
        obs, solverView, discretization);
end

function task = make_observable_pv_guardrail_task()
    sys = make_harmonic_system();
    obs = {@observable_nonlinear_mix};
    params = 1;
    PV = struct('name', 'obs', 'idx', 1);
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    [obs, params, t, x] = canonicalize_trajectory_fixture(obs, params, t, x, 3, PV);
    discretization = fmam_state_defaults.defaultDiscretization();
    extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs, params, t, x, 3, PV, ...
        discretization, extremaSearch);
    derived = fmam_state_ops.buildDerivedView( ...
        obs, solverView, discretization);
    item = struct('prop', 'obsAmp', 'idx', 1, 'target', derived.obsAmp(1));

    task = FMAM_ODE(sys, obs, solverView, item, 1, 0.25, 1e-6, ...
        'derivatives', derivatives);
    task.psiUpdateMode = true;
    task.refreshPsiModeReferences();
    task.needLog = false;
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

function values = current_guardrail_target_values(task)
    solverView = task.exportSolverView();
    derived = task.exportDerivedView();
    values = zeros(1, 5);
    values(1) = solverView.params(2);
    values(2) = solverView.p_Psi(1);
    values(3) = derived.varAmp(1);
    values(4) = derived.obsAmp(1);
    values(5) = derived.varPhase(1, 2);
end

function stat = make_guardrail_state(obs, params)
    PV = struct('name', 'var', 'idx', 1);
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    stat = state(obs, params, t, x, 3, PV);
end

function residual = psi_reference_residual(task)
    errBound = task.errBound;
    task.errBound = inf;
    cleanup = onCleanup(@() restore_err_bound(task, errBound));
    result = task.oneIter();
    M = task.truncationOrder;
    tailLength = max(0,2 * M - 2);
    if tailLength == 0
        residual = zeros(0, 1);
        return
    end
    rowEnd = numel(result.linearResidual) - 1;
    residual = result.linearResidual(rowEnd - tailLength + 1:rowEnd);
end

function restore_err_bound(task, errBound)
    task.errBound = errBound;
end

function verifyPropertyWriteFails(testCase, setter)
    didFail = false;
    try
        setter();
    catch ME
        didFail = true;
        verifyNotEmpty(testCase, ME.identifier);
    end

    verifyTrue(testCase, didFail);
end

function assign_continuation_options(task)
    task.continuationOptions = struct('predictorMode', 'constant');
end

function mutate_solver_view(task, editor)
    view = task.exportSolverView();
    view = editor(view);
    task.loadSolverView(view);
end

function vec = replace_at(vec, idx, value)
    vec(idx) = value;
end

function mat = bump_matrix_entry(mat, rowIdx, colIdx, delta)
    mat(rowIdx, colIdx) = mat(rowIdx, colIdx) + delta;
end

function view = with_field(view, fieldName, value)
    view.(fieldName) = value;
end

function view = apply_psi_delta(view, delta_p, delta_q)
    padded_p = zeros(size(view.p_Psi));
    padded_q = zeros(size(view.q_Psi));
    padded_p(1:min(numel(delta_p), numel(padded_p))) = delta_p(1:min(numel(delta_p), numel(padded_p)));
    padded_q(1:min(numel(delta_q), numel(padded_q))) = delta_q(1:min(numel(delta_q), numel(padded_q)));
    view.p_Psi = view.p_Psi + padded_p;
    view.q_Psi = view.q_Psi + padded_q;
end

function assign_items_perturb(task)
    task.items_perturb = struct('prop', 'params', 'idx', 1, 'target', 2);
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

function z = observable_x1(y)
    z = y(:, 1);
end

function z = observable_nonlinear_mix(y)
    z = y(:, 1) + 0.2 * y(:, 1).^2 + 0.1 * y(:, 2);
end
