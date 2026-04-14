function tests = test_fmam_ode_transaction
%TEST_FMAM_ODE_TRANSACTION Regression coverage for transactional continuation.

    tests = functiontests(localfunctions);
end

function testConstructorRequiresExternalDerivatives(testCase)
    sys = make_harmonic_system();
    obs = {};
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);
    items_per = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1));

    verifyError(testCase, @() FMAM_ODE(sys, obs, solverView, items_per, 1, 1e-6, 1e-6), ...
        'FMAM_ODE:MissingDerivatives');
end

function testFitReturnsStructuredResult(testCase)
    task = make_period_target_task(0, 1e-6, 5);

    result = task.fit();

    verifyTrue(testCase, isstruct(result));
    verifyTrue(testCase, isfield(result, 'converged'));
    verifyTrue(testCase, isfield(result, 'iterations'));
    verifyTrue(testCase, isfield(result, 'finalError'));
    verifyTrue(testCase, isfield(result, 'message'));
    verifyTrue(testCase, isfield(result, 'history'));
    verifyTrue(testCase, isfield(result, 'linearResidualNorm'));
    verifyTrue(testCase, result.converged);
    verifyGreaterThanOrEqual(testCase, result.iterations, 0);
end

function testFitResultContractRemainsStable(testCase)
    task = make_period_target_task(1e-3, 1e-6, 5);

    result = task.fit();
    expectedFields = {'converged', 'iterations', 'finalError', 'message', 'history', ...
        'linearResidualNorm', 'linearResidual', 'stepAccepted', 'scalarError', 'objective'};

    verifyEqual(testCase, sort(fieldnames(result)), sort(expectedFields(:)));
    verifyClass(testCase, result.converged, 'logical');
    verifyTrue(testCase, isscalar(result.iterations) && isnumeric(result.iterations));
    verifyTrue(testCase, isrow(result.finalError));
    verifyEqual(testCase, numel(result.finalError), numel(task.res()));
    verifyTrue(testCase, ischar(result.message) || (isstring(result.message) && isscalar(result.message)));
    verifyTrue(testCase, isscalar(result.linearResidualNorm) && isnumeric(result.linearResidualNorm));
    verifyTrue(testCase, iscolumn(result.linearResidual) && isnumeric(result.linearResidual));
    verifyTrue(testCase, isscalar(result.stepAccepted) && islogical(result.stepAccepted));
    verifyTrue(testCase, isscalar(result.scalarError) && isnumeric(result.scalarError));
    verifyTrue(testCase, isscalar(result.objective) && isnumeric(result.objective));
end

function testFitAndOneIterExposeGenericSolverDiagnostics(testCase)
    task = make_period_target_task(1e-3, 1e-6, 5);

    fitResult = task.fit();

    verifyGreaterThanOrEqual(testCase, numel(fitResult.history), 1);
    verifyHistoryShape(testCase, fitResult.history);

    task = make_period_target_task(1e-3, 1e-6, 5);
    iterResult = task.oneIter();

    verifyTrue(testCase, isstruct(iterResult));
    verifyTrue(testCase, isfield(iterResult, 'history'));
    verifyTrue(testCase, isfield(iterResult, 'stepAccepted'));
    verifyLessThanOrEqual(testCase, numel(iterResult.history), 1);
    verifyEqual(testCase, iterResult.iterations, numel(iterResult.history));
    if ~isempty(iterResult.history)
        verifyHistoryShape(testCase, iterResult.history);
    end
    verifyTrue(testCase, iterResult.stepAccepted);
end

function testOneIterBacktracksWhenPsiValidatorRejectsFullStep(testCase)
    task = make_period_target_task(-1.2, 1e-6, 1, ...
        'continuationOptions', struct('initialLambdaStep', 1));
    task.autostepsize();
    task.perturb();

    result = task.oneIter();
    finalView = task.exportSolverView();
    phi = (0:state.LphiConst-1)' * 2 * pi / state.LphiConst;
    psiMin = min(fmam_state_ops.evaluateTrigSeries(phi, finalView.p_Psi, finalView.q_Psi));

    verifyEqual(testCase, result.iterations, 1);
    verifyEqual(testCase, numel(result.history), 1);
    verifyTrue(testCase, result.stepAccepted);
    verifyFalse(testCase, result.converged);
    verifyLessThan(testCase, result.history.acceptedScale, 1);
    verifyGreaterThan(testCase, result.history.backtracks, 0);
    verifyGreaterThan(testCase, psiMin, 0);
end

function testFailedStepRollsBackAcceptedState(testCase)
    delta = 1e-3;
    task = make_period_target_task(delta, delta, 0);
    acceptedSnapshot = task.exportSolverView();
    acceptedTarget = task.items_per_curr;

    task.step();

    verifyEqual(testCase, task.items_per_curr, acceptedTarget, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.logs, []);
    verifySolverViewSnapshot(testCase, task.exportSolverView(), acceptedSnapshot);
end

function testSuccessfulStepCommitsTargetAndLogs(testCase)
    delta = 1e-3;
    task = make_period_target_task(delta, delta / 10, 20, ...
        'continuationOptions', struct('initialSteps', 1));
    acceptedTarget = task.items_per_curr;
    targetValue = task.items_perturb.target;

    task.step();

    verifyGreaterThan(testCase, task.items_per_curr, acceptedTarget);
    verifyEqual(testCase, task.items_per_curr, targetValue, 'AbsTol', 1e-12);
    verifyEqual(testCase, numel(task.logs), 1);
    verifyEqual(testCase, task.logs.lambda, 1, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.logs.predictorUsed, 'constant');
    finalView = task.exportSolverView();
    verifyEqual(testCase, finalView.p_Psi(1), targetValue, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.logs.p_Psi_idx_1, targetValue, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.logs.accepted_p_Psi_idx_1, finalView.p_Psi(1), 'AbsTol', 1e-12);
end

function testConditionStopDiscardRestoresAcceptedState(testCase)
    delta = 1e-3;
    task = make_mock_period_task(delta, inf, struct( ...
        'initialSteps', 1, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant', ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', 2));
    acceptedSnapshot = task.exportSolverView();
    acceptedTarget = task.items_per_curr;

    task.step();

    verifyEqual(testCase, task.items_per_curr, acceptedTarget, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.logs, []);
    verifySolverViewSnapshot(testCase, task.exportSolverView(), acceptedSnapshot);
    verifyFalse(testCase, task.continuationStatus.completed);
    verifyEqual(testCase, task.continuationStatus.lambda, 0, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.continuationStatus.reason, 'condition_stop');
    verifyEqual(testCase, task.continuationStatus.triggerValue, 1, 'AbsTol', 1e-12);
end

function testRejectedTrialShrinksAndRetries(testCase)
    task = make_mock_period_task(0.4, 0.25, struct( ...
        'initialSteps', 1, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant'));

    task.step();

    verifyEqual(testCase, numel(task.logs), 2);
    verifyEqual(testCase, [task.logs.lambda], [0.5, 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, task.items_per_curr, task.items_perturb.target, 'AbsTol', 1e-12);
end

function testRejectedRetryMatchesCleanAcceptedState(testCase)
    retryTask = make_mock_period_task(0.4, 0.25, struct( ...
        'initialSteps', 1, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant'));
    baselineTask = make_mock_period_task(0.4, inf, struct( ...
        'initialSteps', 2, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant'));

    retryTask.step();
    baselineTask.step();

    verifyEqual(testCase, [retryTask.logs.lambda], [0.5, 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, [baselineTask.logs.lambda], [0.5, 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, retryTask.items_per_curr, baselineTask.items_per_curr, 'AbsTol', 1e-12);
    verifySolverViewSnapshot(testCase, retryTask.exportSolverView(), baselineTask.exportSolverView());
end

function testConvergedStateWithInvalidPsiRollsBackAcceptedState(testCase)
    delta = 1e-3;
    task = make_mock_period_task(delta, inf, struct( ...
        'initialSteps', 1, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'maxFailures', 0));
    task.mockInvalidPsiAboveStep = 0;
    acceptedSnapshot = task.exportSolverView();
    acceptedTarget = task.items_per_curr;

    task.step();

    verifyEqual(testCase, task.items_per_curr, acceptedTarget, 'AbsTol', 1e-12);
    verifyEqual(testCase, task.logs, []);
    verifySolverViewSnapshot(testCase, task.exportSolverView(), acceptedSnapshot);
end

function testInvalidPsiConvergedStateShrinksAndRetries(testCase)
    task = make_mock_period_task(0.4, inf, struct( ...
        'initialSteps', 1, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant'));
    task.mockInvalidPsiAboveStep = 0.25;

    task.step();

    verifyEqual(testCase, numel(task.logs), 2);
    verifyEqual(testCase, [task.logs.lambda], [0.5, 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, task.items_per_curr, task.items_perturb.target, 'AbsTol', 1e-12);
end

function testAutoPredictorUsesHermiteAfterTwoAcceptedPoints(testCase)
    task = make_mock_period_task(0.3, inf, struct( ...
        'initialSteps', 2, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1));

    task.step();

    verifyEqual(testCase, numel(task.logs), 2);
    verifyTrue(testCase, any(strcmp(task.logs(2).predictorUsed, {'hermite', 'secant'})));
end

function testAutoPredictorUsesQuadraticAfterThreeAcceptedPoints(testCase)
    task = make_mock_period_task(0.3, inf, struct( ...
        'initialSteps', 3, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1));

    task.step();

    verifyEqual(testCase, numel(task.logs), 3);
    verifyEqual(testCase, task.logs(3).predictorUsed, 'quadratic');
    verifyEqual(testCase, task.logs(3).predictorFallbackCount, 0);
end

function testPredictorValidationFallbackDoesNotPolluteAcceptedState(testCase)
    fallbackTask = make_mock_period_task(0.3, inf, struct( ...
        'initialSteps', 3, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'quadratic', ...
        'predictorStepGrowthLimit', 1e-12));
    baselineTask = make_mock_period_task(0.3, inf, struct( ...
        'initialSteps', 3, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant'));

    fallbackTask.step();
    baselineTask.step();

    verifyEqual(testCase, numel(fallbackTask.logs), 3);
    verifyGreaterThan(testCase, max([fallbackTask.logs.predictorFallbackCount]), 0);
    verifyEqual(testCase, fallbackTask.items_per_curr, baselineTask.items_per_curr, 'AbsTol', 1e-12);
    verifySolverViewSnapshot(testCase, fallbackTask.exportSolverView(), baselineTask.exportSolverView());
end

function testHigherOrderPredictorRetryRestoresCurrentAcceptedState(testCase)
    retryTask = make_mock_period_task(1, 0.3, struct( ...
        'initialSteps', 8, ...
        'growthFactor', 2, ...
        'shrinkFactor', 1, ...
        'backtrackingFactor', 0.5, ...
        'predictorMode', 'quadratic', ...
        'predictorStepGrowthLimit', 1e-12));
    baselineTask = make_mock_period_task(1, 0.3, struct( ...
        'initialSteps', 8, ...
        'growthFactor', 2, ...
        'shrinkFactor', 1, ...
        'backtrackingFactor', 0.5, ...
        'predictorMode', 'constant'));

    retryTask.step();
    baselineTask.step();

    verifyGreaterThanOrEqual(testCase, numel(retryTask.logs), 4);
    verifyGreaterThan(testCase, max([retryTask.logs.predictorFallbackCount]), 0);
    verifyEqual(testCase, [retryTask.logs.lambda], [baselineTask.logs.lambda], 'AbsTol', 1e-12);
    verifyEqual(testCase, retryTask.items_per_curr, baselineTask.items_per_curr, 'AbsTol', 1e-12);
    verifySolverViewSnapshot(testCase, retryTask.exportSolverView(), baselineTask.exportSolverView());
end

function testSuccessfulContinuationRestoresLastAcceptedSolverState(testCase)
    task = make_mock_period_task(0.4, inf, struct( ...
        'initialSteps', 2, ...
        'growthFactor', 1, ...
        'shrinkFactor', 1, ...
        'predictorMode', 'constant'));
    targetValue = task.items_perturb.target;

    task.step();

    verifyEqual(testCase, [task.logs.lambda], [0.5, 1], 'AbsTol', 1e-12);
    verifyTrue(testCase, task.continuationStatus.completed);
    verifyEqual(testCase, task.continuationStatus.reason, 'target_reached');
    verifyEqual(testCase, task.items_per_curr, targetValue, 'AbsTol', 1e-12);
    finalView = task.exportSolverView();
    verifyEqual(testCase, finalView.p_Psi(1), targetValue, 'AbsTol', 1e-12);
    verifyEqual(testCase, finalView.params(1), 1 / targetValue, 'AbsTol', 1e-12);
end

function testWrappedPhaseTargetUsesShortestEquivalentPath(testCase)
    task = make_mock_phase_task(2, 2 * pi - 0.1, struct('initialSteps', 1));
    solverView = task.exportSolverView();
    initialPhi = solverView.varPhiMax(2);

    task.step();

    verifyEqual(testCase, numel(task.logs), 1);
    verifyEqual(testCase, task.items_per_curr, initialPhi - 0.1, 'AbsTol', 1e-12);
    verifyLessThan(testCase, task.items_per_curr, task.items_perturb.target);
end

function task = make_period_target_task(delta, accuracy, maxIterations, varargin)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);

    items_per = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1) + delta);
    maxStep = max(abs(delta), 1e-6);
    newtonOptions = struct('maxIterations', maxIterations);
    task = FMAM_ODE(sys, obs, solverView, items_per, 1, maxStep, accuracy, ...
        'derivatives', derivatives, 'newtonOptions', newtonOptions, varargin{:});
    task.needLog = true;
end

function task = make_mock_period_task(delta, failAboveStep, continuationOptions)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);

    items_per = struct('prop', 'p_Psi', 'idx', 1, 'target', solverView.p_Psi(1) + delta);
    maxStep = max(abs(delta), 1e-6);
    task = MockFMAMContinuationTask(sys, obs, solverView, items_per, 1, maxStep, 1e-6, ...
        'derivatives', derivatives, ...
        'continuationOptions', continuationOptions);
    task.needLog = true;
    task.mockFailAboveStep = failAboveStep;
end

function task = make_mock_phase_task(varIdx, rawShift, continuationOptions)
    sys = make_harmonic_system();
    obs = {};
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    PV = struct('name', 'var', 'idx', 1);
    solverView = make_reference_solver_view(obs, 1, 3, PV);

    items_per = struct('prop', 'varPhiMax', 'idx', varIdx, ...
        'target', solverView.varPhiMax(varIdx) + rawShift);
    task = MockFMAMContinuationTask(sys, obs, solverView, items_per, 1, abs(rawShift), 1e-6, ...
        'derivatives', derivatives, ...
        'continuationOptions', continuationOptions);
    task.needLog = true;
end

function solverView = make_reference_solver_view(obs, params, M, PV)
    [t, x] = reference_trajectory();
    [obs, params, t, x] = canonicalize_trajectory_fixture(obs, params, t, x, M, PV);
    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs, params, t, x, M, PV, ...
        state.Lconst, state.LphiConst, ...
        state.countMax, state.errMax);
end

function [t, x] = reference_trajectory()
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
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

function verifySolverViewSnapshot(testCase, view, snapshot)
    fields = {'params', 'p_Psi', 'q_Psi', 'p_var', 'q_var', ...
        'varPhiMax', 'varPhiMin', 'obsPhiMax', 'obsPhiMin'};

    for i = 1:numel(fields)
        fieldName = fields{i};
        verifyEqual(testCase, view.(fieldName), snapshot.(fieldName), 'AbsTol', 1e-12);
    end
end

function verifyHistoryShape(testCase, history)
    requiredFields = {'iteration', 'objective', 'residualNorm', 'stepNorm', ...
        'acceptedScale', 'lambda', 'backtracks', 'accepted', 'maxError', ...
        'solver', 'conditionEstimate'};

    for i = 1:numel(requiredFields)
        verifyTrue(testCase, isfield(history, requiredFields{i}));
    end
end
