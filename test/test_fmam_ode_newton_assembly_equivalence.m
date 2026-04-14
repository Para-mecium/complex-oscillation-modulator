function tests = test_fmam_ode_newton_assembly_equivalence
%TEST_FMAM_ODE_NEWTON_ASSEMBLY_EQUIVALENCE Freeze Newton assembly behavior.

    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(testDir);
    addpath(rootDir);
    testCase.TestData.rootDir = rootDir;
end

function testDirectCoefficientTargetMatchesSharedHelper(testCase)
    task = make_task_fixture(struct('prop', 'p_Psi', 'idx', 1), {}, ...
        struct('name', 'var', 'idx', 1), false);

    verifyAssemblyMatchesSharedHelper(testCase, task);
end

function testVariableExtremaTargetMatchesSharedHelper(testCase)
    task = make_task_fixture(struct('prop', 'varAmp', 'idx', 1), {}, ...
        struct('name', 'var', 'idx', 1), false);

    verifyAssemblyMatchesSharedHelper(testCase, task);
end

function testObservableExtremaTargetMatchesSharedHelper(testCase)
    obs = {@observable_x1};
    task = make_task_fixture(struct('prop', 'obsAmp', 'idx', 1), obs, ...
        struct('name', 'obs', 'idx', 1), true);

    verifyAssemblyMatchesSharedHelper(testCase, task);
end

function testObservablePVGaugeRowsMatchFiniteDifference(testCase)
    obs = {@observable_nonlinear_mix};
    task = make_current_target_task(struct('prop', 'obsAmp', 'idx', 1), obs, ...
        struct('name', 'obs', 'idx', 1), true);
    ctx = make_assembly_context(task);

    [A, ~, indexMap] = assemble_newton_linear_system(ctx);
    gaugeRows = observable_pv_gauge_row_indices(ctx, size(A, 1));

    fdP = finite_difference_observable_gauge_column(ctx, ...
        @(view, epsVal) set_p_var_entry(view, 3, 1, epsVal));
    fdQ = finite_difference_observable_gauge_column(ctx, ...
        @(view, epsVal) set_q_var_entry(view, 1, 1, epsVal));

    verifyEqual(testCase, A(gaugeRows, indexMap.p_var(3, 1)), fdP, ...
        'AbsTol', 1e-6, 'RelTol', 1e-5);
    verifyEqual(testCase, A(gaugeRows, indexMap.q_var(1, 1)), fdQ, ...
        'AbsTol', 1e-6, 'RelTol', 1e-5);
end

function testObservableExtremaRowsMatchFiniteDifference(testCase)
    obs = {@observable_nonlinear_mix};
    task = make_current_target_task(struct('prop', 'obsAmp', 'idx', 1), obs, ...
        struct('name', 'obs', 'idx', 1), true);
    ctx = make_assembly_context(task);

    [A, res, indexMap] = assemble_newton_linear_system(ctx);
    rowMax = observable_extrema_row_index(ctx, 1, 'max');
    rowMin = observable_extrema_row_index(ctx, 1, 'min');

    verifyObservableExtremaRowFiniteDifference(testCase, ctx, A, res, indexMap, rowMax, 1, 'max');
    verifyObservableExtremaRowFiniteDifference(testCase, ctx, A, res, indexMap, rowMin, 1, 'min');
end

function testVarPhaseTargetMatchesSharedHelper(testCase)
    task = make_task_fixture(struct('prop', 'varPhase', 'idx', [1, 2]), {}, ...
        struct('name', 'var', 'idx', 1), false);

    verifyAssemblyMatchesSharedHelper(testCase, task);
end

function testOneIterMatchesAssemblyResidualForFrozenPsiMode(testCase)
    task = make_current_target_task(struct('prop', 'p_Psi', 'idx', 1), {}, ...
        struct('name', 'var', 'idx', 1), false);

    verifyOneIterMatchesAssemblyResidual(testCase, task);
end

function testOneIterMatchesAssemblyResidualForUpdatedPsiMode(testCase)
    task = make_current_target_task(struct('prop', 'p_Psi', 'idx', 1), {}, ...
        struct('name', 'var', 'idx', 1), true);

    verifyOneIterMatchesAssemblyResidual(testCase, task);
end

function verifyAssemblyMatchesSharedHelper(testCase, task)
    ctx = make_assembly_context(task);

    [A, res, indexMap, unknowns] = assemble_newton_linear_system(ctx);

    verifyEqual(testCase, unknowns, {'params', 'p_Psi', 'q_Psi', 'p_var', 'q_var', ...
        'varPhiMax', 'varPhiMin', 'obsPhiMax', 'obsPhiMin'});
    verifyTrue(testCase, all(isfield(indexMap, unknowns)));
    verifyEqual(testCase, size(A,1), numel(res));
    maxIndex = max(cellfun(@(idx) max([idx(:); 0]), struct2cell(indexMap)));
    verifyEqual(testCase, size(A,2), maxIndex);
    verifyFalse(testCase, any(isnan(A(:))));
    verifyFalse(testCase, any(isnan(res(:))));

    verifyTargetRowsMatchSharedHelper(testCase, ctx, A, res, indexMap);
    verifyExtremaRowsRespectSharedHelper(testCase, ctx, A, res, indexMap);
end

function verifyOneIterMatchesAssemblyResidual(testCase, task)
    acceptedTarget = task.items_per_curr;

    result = task.oneIter();
    ctx = make_assembly_context(task);
    [~, expectedResidual] = assemble_newton_linear_system(ctx);
    expectedErrorVec = reshape(task.res(), 1, []);

    verifyEqual(testCase, result.iterations, numel(result.history));
    verifyLessThanOrEqual(testCase, numel(result.history), 1);
    if ~isempty(result.history)
        verifyNewtonHistoryShape(testCase, result.history);
    end
    verifyEqual(testCase, result.linearResidual(:), expectedResidual(:), 'AbsTol', 1e-10);
    verifyEqual(testCase, result.linearResidualNorm, norm(expectedResidual, 2), 'AbsTol', 1e-10);
    verifyEqual(testCase, reshape(result.errorVec, 1, []), expectedErrorVec, 'AbsTol', 1e-10);
    verifyEqual(testCase, result.scalarError, max(expectedErrorVec), 'AbsTol', 1e-10);
    verifyEqual(testCase, task.items_per_curr, acceptedTarget, 'AbsTol', 1e-12);
end

function task = make_task_fixture(itemTemplate, obs, PV, isPsiUpdated)
    sys = make_harmonic_system();
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    [solverView, derived] = make_reference_views(obs, 1, 3, PV);

    item = itemTemplate;
    targetCtx = make_target_rule_context(solverView, derived, obs, numel(sys));
    values = make_target_value_struct(solverView);
    item.target = fmam_target_rules('current_value', targetCtx, item, values);
    delta = 1e-3;
    item.target = item.target + delta;

    task = FMAM_ODE(sys, obs, solverView, item, 1, delta, 1e-6, ...
        'derivatives', derivatives);
    task.isPsiUpdated = isPsiUpdated;
    task.autostepsize();
    task.perturb();
end

function task = make_current_target_task(itemTemplate, obs, PV, isPsiUpdated)
    sys = make_harmonic_system();
    derivatives = build_symbolic_derivatives(sys, obs, 1);
    [solverView, derived] = make_reference_views(obs, 1, 3, PV);

    item = itemTemplate;
    targetCtx = make_target_rule_context(solverView, derived, obs, numel(sys));
    values = make_target_value_struct(solverView);
    item.target = fmam_target_rules('current_value', targetCtx, item, values);

    task = FMAM_ODE(sys, obs, solverView, item, 1, 1e-3, 1e-6, ...
        'derivatives', derivatives);
    task.isPsiUpdated = isPsiUpdated;
end

function verifyTargetRowsMatchSharedHelper(testCase, ctx, A, res, indexMap)
    firstTargetRow = first_target_row_index(ctx);
    assembly = make_target_row_assembly(ctx, indexMap, size(A, 2));
    for j = 1:numel(ctx.items_perturb)
        [targetRow, targetResidual] = fmam_target_rules('target_row', ...
            ctx.targetRuleCtx, ctx.items_perturb(j), ctx.targetCurr(j), assembly);
        rowIdx = firstTargetRow + j - 1;
        verifyEqual(testCase, A(rowIdx, :), targetRow);
        verifyEqual(testCase, res(rowIdx), targetResidual);
    end
end

function verifyExtremaRowsRespectSharedHelper(testCase, ctx, A, res, indexMap)
    rowIdx = first_extrema_row_index(ctx);

    for i = 1:ctx.dimVar
        [needMax, needMin] = fmam_target_rules('needs_variable_extrema', ctx.targetRuleCtx, ctx.items_perturb, i);
        verifyExtremaGuardRow(testCase, A(rowIdx, :), res(rowIdx), ...
            indexMap.varPhiMax(i), needMax);
        rowIdx = rowIdx + 1;

        verifyExtremaGuardRow(testCase, A(rowIdx, :), res(rowIdx), ...
            indexMap.varPhiMin(i), needMin);
        rowIdx = rowIdx + 1;
    end

    for i = 1:ctx.dimObs
        [needMax, needMin] = fmam_target_rules('needs_observable_extrema', ctx.targetRuleCtx, ctx.items_perturb, i);
        verifyExtremaGuardRow(testCase, A(rowIdx, :), res(rowIdx), ...
            indexMap.obsPhiMax(i), needMax);
        rowIdx = rowIdx + 1;

        verifyExtremaGuardRow(testCase, A(rowIdx, :), res(rowIdx), ...
            indexMap.obsPhiMin(i), needMin);
        rowIdx = rowIdx + 1;
    end
end

function verifyExtremaGuardRow(testCase, row, residual, phiIdx, isActive)
    identityRow = zeros(size(row));
    identityRow(phiIdx) = 1;

    if isActive
        verifyNotEqual(testCase, row, identityRow);
    else
        verifyEqual(testCase, row, identityRow);
        verifyEqual(testCase, residual, 0);
    end
end

function rowIdx = first_extrema_row_index(ctx)
    rowIdx = ctx.dimVar * (2 * ctx.truncationOrder + 1);
    if ~ctx.isPsiUpdated
        rowIdx = rowIdx + 1;
    end
    rowIdx = rowIdx + 1;
end

function rowIdx = first_target_row_index(ctx)
    rowIdx = first_extrema_row_index(ctx);
    rowIdx = rowIdx + 2 * ctx.dimVar + 2 * ctx.dimObs;
    rowIdx = rowIdx + (ctx.dimParams - numel(ctx.items_controlled));
end

function ctx = make_target_rule_context(solverView, derived, obs, dimVar)
    ctx = struct();
    ctx.obs = obs;
    ctx.dimVar = dimVar;
    ctx.dimObs = numel(obs);
    ctx.dimParams = numel(solverView.params);
    ctx.propertySizes = struct( ...
        'p_Psi', size(solverView.p_Psi), ...
        'q_Psi', size(solverView.q_Psi), ...
        'p_var', size(solverView.p_var), ...
        'q_var', size(solverView.q_var), ...
        'varPhiMax', size(solverView.varPhiMax), ...
        'varPhiMin', size(solverView.varPhiMin), ...
        'obsPhiMax', size(solverView.obsPhiMax), ...
        'obsPhiMin', size(solverView.obsPhiMin), ...
        'params', size(solverView.params));
    ctx.solver = solverView;
    ctx.derived = derived;
end

function values = make_target_value_struct(solverView)
    values = struct();
    values.parameters = solverView.params;
    values.p_Psi = solverView.p_Psi;
    values.q_Psi = solverView.q_Psi;
    values.p_var = solverView.p_var;
    values.q_var = solverView.q_var;
end

function [solverView, derived] = make_reference_views(obs, omega, M, PV)
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    [obs, omega, t, x] = canonicalize_trajectory_fixture(obs, omega, t, x, M, PV);
    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs, omega, t, x, M, PV, ...
        state.Lconst, state.LphiConst, ...
        state.countMax, state.errMax);
    derived = fmam_state_ops.buildDerivedView( ...
        obs, solverView, state.LphiConst, state.Lconst, ...
        state.countMax, state.errMax);
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

function assembly = make_target_row_assembly(ctx, indexMap, numUnknown)
    assembly = struct();
    assembly.indexMap = indexMap;
    assembly.parameters = ctx.params;
    assembly.p_Psi = ctx.p_Psi;
    assembly.q_Psi = ctx.q_Psi;
    assembly.p_var = ctx.p_var;
    assembly.q_var = ctx.q_var;
    assembly.Derivative_obs = ctx.derivatives.obs;
    assembly.numUnknown = numUnknown;
end

function ctx = make_assembly_context(task)
    solverView = task.exportSolverView();
    derived = task.exportDerivedView();
    ctx = struct();
    ctx.sys = task.sys;
    ctx.obs = task.obs;
    ctx.derivatives = task.derivatives;
    ctx.L = task.Lconst;
    ctx.dimVar = size(solverView.p_var, 2);
    ctx.dimObs = numel(task.obs);
    ctx.truncationOrder = size(solverView.q_var, 1);
    ctx.dimParams = numel(solverView.params);
    ctx.params = solverView.params;
    ctx.p_Psi = solverView.p_Psi;
    ctx.q_Psi = solverView.q_Psi;
    ctx.p_var = solverView.p_var;
    ctx.q_var = solverView.q_var;
    ctx.varPhiMax = solverView.varPhiMax;
    ctx.varPhiMin = solverView.varPhiMin;
    ctx.obsPhiMax = solverView.obsPhiMax;
    ctx.obsPhiMin = solverView.obsPhiMin;
    ctx.PV = solverView.PV;
    ctx.items_perturb = task.items_perturb;
    ctx.items_controlled = task.items_controlled;
    ctx.targetCurr = task.items_per_curr;
    ctx.isPsiUpdated = task.isPsiUpdated;
    ctx.p_Psi_init = solverView.p_Psi;
    ctx.q_Psi_init = solverView.q_Psi;
    ctx.targetRuleCtx = make_target_rule_context(solverView, derived, task.obs, ctx.dimVar);
end

function verifyNewtonHistoryShape(testCase, history)
    requiredFields = {'iteration', 'objective', 'residualNorm', 'stepNorm', ...
        'acceptedScale', 'lambda', 'backtracks', 'accepted', 'maxError', ...
        'solver', 'conditionEstimate', 'directConditionEstimate'};

    for i = 1:numel(requiredFields)
        verifyTrue(testCase, isfield(history, requiredFields{i}));
    end
end

function verifyObservableExtremaRowFiniteDifference(testCase, ctx, A, res, indexMap, rowIdx, obsIdx, kind)
    M = ctx.truncationOrder;
    epsVal = 1e-7;
    solverView = solver_view_from_context(ctx);
    verifyEqual(testCase, res(rowIdx), ...
        observable_extrema_row_residual(ctx, solverView, obsIdx, kind), ...
        'AbsTol', 1e-12);

    for varIdx = 1:ctx.dimVar
        for pIdx = 1:(M + 1)
            colIdx = indexMap.p_var(pIdx, varIdx);
            fd = finite_difference_observable_extrema_row(ctx, obsIdx, kind, ...
                @(view) perturb_p_var(view, pIdx, varIdx, epsVal), epsVal);
            verifyEqual(testCase, -A(rowIdx, colIdx), fd, 'AbsTol', 1e-6, 'RelTol', 1e-5);
        end

        for qIdx = 1:M
            colIdx = indexMap.q_var(qIdx, varIdx);
            fd = finite_difference_observable_extrema_row(ctx, obsIdx, kind, ...
                @(view) perturb_q_var(view, qIdx, varIdx, epsVal), epsVal);
            verifyEqual(testCase, -A(rowIdx, colIdx), fd, 'AbsTol', 1e-6, 'RelTol', 1e-5);
        end
    end

    if strcmp(kind, 'max')
        colIdx = indexMap.obsPhiMax(obsIdx);
        fd = finite_difference_observable_extrema_row(ctx, obsIdx, kind, ...
            @(view) perturb_obs_phi(view, obsIdx, epsVal, kind), epsVal);
    else
        colIdx = indexMap.obsPhiMin(obsIdx);
        fd = finite_difference_observable_extrema_row(ctx, obsIdx, kind, ...
            @(view) perturb_obs_phi(view, obsIdx, epsVal, kind), epsVal);
    end
    verifyEqual(testCase, -A(rowIdx, colIdx), fd, 'AbsTol', 1e-6, 'RelTol', 1e-5);
end

function rowIdx = observable_pv_gauge_row_indices(ctx, numRows)
    M = ctx.truncationOrder;
    rowIdx = (numRows - (2 * M - 1) + 1):numRows;
end

function rowIdx = observable_extrema_row_index(ctx, obsIdx, kind)
    rowIdx = first_extrema_row_index(ctx) + 2 * ctx.dimVar + 2 * (obsIdx - 1);
    if strcmp(kind, 'min')
        rowIdx = rowIdx + 1;
    end
end

function fd = finite_difference_observable_gauge_column(ctx, editor)
    solverView = solver_view_from_context(ctx);
    gauge0 = observable_primary_gauge(solverView, ctx.obs, ctx.PV, ctx.L);
    epsVal = 1e-7;
    solverViewPerturbed = editor(solverView, epsVal);
    gauge1 = observable_primary_gauge(solverViewPerturbed, ctx.obs, ctx.PV, ctx.L);
    fd = (gauge1 - gauge0) / epsVal;
end

function fd = finite_difference_observable_extrema_row(ctx, obsIdx, kind, editor, epsVal)
    solverView = solver_view_from_context(ctx);
    row0 = observable_extrema_row_residual(ctx, solverView, obsIdx, kind);
    solverViewPerturbed = editor(solverView);
    row1 = observable_extrema_row_residual(ctx, solverViewPerturbed, obsIdx, kind);
    fd = (row1 - row0) / epsVal;
end

function residual = observable_extrema_row_residual(ctx, solverView, obsIdx, kind)
    if strcmp(kind, 'max')
        phi = solverView.obsPhiMax(obsIdx);
    else
        phi = solverView.obsPhiMin(obsIdx);
    end
    residual = -FMAM_ODE.residue_phi_obs( ...
        ctx.derivatives.obs, solverView.p_var, solverView.q_var, phi, obsIdx);
end

function view = solver_view_from_context(ctx)
    view = struct( ...
        'params', ctx.params, ...
        'p_Psi', ctx.p_Psi, ...
        'q_Psi', ctx.q_Psi, ...
        'p_var', ctx.p_var, ...
        'q_var', ctx.q_var, ...
        'varPhiMax', ctx.varPhiMax, ...
        'varPhiMin', ctx.varPhiMin, ...
        'obsPhiMax', ctx.obsPhiMax, ...
        'obsPhiMin', ctx.obsPhiMin, ...
        'PV', ctx.PV);
end

function gauge = observable_primary_gauge(view, obs, PV, L)
    M = size(view.q_var, 1);
    phi = (0:L-1).' * 2 * pi / L;
    [vc, vs] = FMAM_ODE.Vec_CS(phi, M, L);
    TS_var = vc * view.p_var + vs * view.q_var;
    coeff = fft(obs{PV.idx}(TS_var)) / L;
    gauge = [real(coeff(3:M+1)); imag(coeff(2:M+1))];
end

function view = set_p_var_entry(view, rowIdx, colIdx, delta)
    view.p_var(rowIdx, colIdx) = view.p_var(rowIdx, colIdx) + delta;
end

function view = set_q_var_entry(view, rowIdx, colIdx, delta)
    view.q_var(rowIdx, colIdx) = view.q_var(rowIdx, colIdx) + delta;
end

function view = perturb_p_var(view, rowIdx, colIdx, delta)
    view.p_var(rowIdx, colIdx) = view.p_var(rowIdx, colIdx) + delta;
end

function view = perturb_q_var(view, rowIdx, colIdx, delta)
    view.q_var(rowIdx, colIdx) = view.q_var(rowIdx, colIdx) + delta;
end

function view = perturb_obs_phi(view, obsIdx, delta, kind)
    if strcmp(kind, 'max')
        view.obsPhiMax(obsIdx) = view.obsPhiMax(obsIdx) + delta;
    else
        view.obsPhiMin(obsIdx) = view.obsPhiMin(obsIdx) + delta;
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

function y = observable_x1(x)
    y = x(:, 1);
end

function y = observable_nonlinear_mix(x)
    y = x(:, 1) + 0.2 * x(:, 1).^2 + 0.1 * x(:, 2);
end
