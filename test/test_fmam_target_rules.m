function tests = test_fmam_target_rules
%TEST_FMAM_TARGET_RULES Direct coverage for shared modulation target semantics.

    tests = functiontests(localfunctions);
end

function testCanonicalizePropertySuccessAndFailure(testCase)
    verifyEqual(testCase, fmam_target_canonicalize("p_psi"), 'p_Psi');
    verifyEqual(testCase, fmam_target_canonicalize('varphase'), 'varPhase');
    verifyCallFails(testCase, @() fmam_target_canonicalize('obsPhase'));
    verifyCallFails(testCase, @() fmam_target_canonicalize(17));
end

function testValidateItemAcceptsSupportedIndicesAndRejectsInvalidOnes(testCase)
    ctx = make_target_ctx({@observable_x1});

    item = struct('prop', 'p_var', 'idx', [2 1], 'target', 0);
    validated = fmam_target_rules('validate_item', ctx, item);
    verifyEqual(testCase, validated.prop, 'p_var');

    verifyCallFails(testCase, @() fmam_target_rules('validate_item', ctx, ...
        struct('prop', 'obsAmp', 'idx', 2, 'target', 0)));
    verifyCallFails(testCase, @() fmam_target_rules('validate_item', ctx, ...
        struct('prop', 'params', 'idx', [1 2], 'target', 0)));
end

function testExtremaNeedsReflectTargetInventory(testCase)
    ctx = make_target_ctx({@observable_x1});
    items(1) = struct('prop', 'varAmp', 'idx', 1, 'target', 0);
    items(2) = struct('prop', 'varMin', 'idx', 1, 'target', 0);
    items(3) = struct('prop', 'obsAmp', 'idx', 1, 'target', 0);

    [needVarMax,needVarMin] = fmam_target_rules('needs_variable_extrema', ctx, items, 1);
    [needObsMax,needObsMin] = fmam_target_rules('needs_observable_extrema', ctx, items, 1);

    verifyTrue(testCase, needVarMax);
    verifyTrue(testCase, needVarMin);
    verifyTrue(testCase, needObsMax);
    verifyTrue(testCase, needObsMin);
end

function testCurrentValueAndLinearIndexUseSharedSemantics(testCase)
    ctx = make_target_ctx({@observable_x1});
    values = make_value_payload(ctx.solver);

    ampItem = struct('prop', 'varAmp', 'idx', 1, 'target', 0);
    phaseItem = struct('prop', 'varPhase', 'idx', [1 2], 'target', 0);
    obsItem = struct('prop', 'obsMax', 'idx', 1, 'target', 0);

    verifyEqual(testCase, fmam_target_rules('current_value', ctx, ampItem, values), ...
        ctx.derived.varAmp(1), ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, fmam_target_rules('current_value', ctx, phaseItem, values), ...
        ctx.derived.varPhase(1,2), ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, fmam_target_rules('current_value', ctx, obsItem, values), ...
        ctx.derived.obsMax(1), ...
        'AbsTol', 1e-12);

    verifyEqual(testCase, fmam_target_rules('linear_index', ctx, 'p_var', [2 1]), 2);
    verifyEqual(testCase, fmam_target_rules('linear_index', ctx, 'params', 3), 3);
end

function testContinuationDeltaWrapsPeriodicTargets(testCase)
    ctx = make_target_ctx({@observable_x1});

    phiItem = struct('prop', 'varPhiMax', 'idx', 1, 'target', 0);
    phaseItem = struct('prop', 'varPhase', 'idx', [1 2], 'target', 0);

    [deltaPhi,isWrappedPhi,periodPhi] = fmam_target_rules('continuation_delta', ctx, ...
        phiItem, 0.1, 2 * pi);
    [deltaPhase,isWrappedPhase,periodPhase] = fmam_target_rules('continuation_delta', ctx, ...
        phaseItem, 0.2, ctx.derived.period + 0.1);

    verifyEqual(testCase, deltaPhi, -0.1, 'AbsTol', 1e-12);
    verifyTrue(testCase, isWrappedPhi);
    verifyEqual(testCase, periodPhi, 2 * pi, 'AbsTol', 1e-12);

    verifyEqual(testCase, deltaPhase, -0.1, 'AbsTol', 1e-12);
    verifyTrue(testCase, isWrappedPhase);
    verifyEqual(testCase, periodPhase, ctx.derived.period, 'AbsTol', 1e-12);
end

function testTargetRowBuildsAmplitudeEquation(testCase)
    ctx = make_target_ctx({@observable_x1});
    assembly = make_var_amp_assembly(ctx.solver);
    item = struct('prop', 'varAmp', 'idx', 1, 'target', 0);
    targetValue = 1.25;

    [row,residual] = fmam_target_rules('target_row', ctx, item, targetValue, assembly);

    [coePMax,coeQMax,coePhiMax] = FMAM_ODE.delta_coe_var_target(assembly.p_var(:,1),assembly.q_var(:,1),ctx.solver.varPhiMax(1));
    [coePMin,coeQMin,coePhiMin] = FMAM_ODE.delta_coe_var_target(assembly.p_var(:,1),assembly.q_var(:,1),ctx.solver.varPhiMin(1));
    expectedRow = zeros(1,assembly.numUnknown);
    expectedRow(assembly.indexMap.p_var(:,1)) = 0.5 * (coePMax - coePMin);
    expectedRow(assembly.indexMap.q_var(:,1)) = 0.5 * (coeQMax - coeQMin);
    expectedRow(assembly.indexMap.varPhiMax(1)) = 0.5 * coePhiMax;
    expectedRow(assembly.indexMap.varPhiMin(1)) = -0.5 * coePhiMin;
    expectedResidual = targetValue - fmam_target_rules('current_value', ctx, item, make_value_payload(ctx.solver));

    verifyEqual(testCase, row, expectedRow, 'AbsTol', 1e-12);
    verifyEqual(testCase, residual, expectedResidual, 'AbsTol', 1e-12);
end

function ctx = make_target_ctx(obs)
    params = ones(1,5);
    stat = make_guardrail_state(obs, params);
    solver = struct( ...
        'params', stat.params, ...
        'p_Psi', stat.p_Psi, ...
        'q_Psi', stat.q_Psi, ...
        'p_var', stat.p_var, ...
        'q_var', stat.q_var, ...
        'varPhiMax', stat.varPhiMax, ...
        'varPhiMin', stat.varPhiMin, ...
        'obsPhiMax', stat.obsPhiMax, ...
        'obsPhiMin', stat.obsPhiMin, ...
        'PV', stat.PV);
    discretization = stat.discretization;
    derived = fmam_state_ops.buildDerivedView( ...
        reshape(obs,1,[]),solver,discretization);
    ctx = struct( ...
        'obs', {reshape(obs,1,[])}, ...
        'solver', solver, ...
        'derived', derived, ...
        'dimVar', 2, ...
        'dimObs', numel(obs), ...
        'dimParams', numel(params), ...
        'propertySizes', struct( ...
            'params', size(stat.params), ...
            'p_Psi', size(stat.p_Psi), ...
            'q_Psi', size(stat.q_Psi), ...
            'p_var', size(stat.p_var), ...
            'q_var', size(stat.q_var), ...
            'varPhiMax', size(stat.varPhiMax), ...
            'varPhiMin', size(stat.varPhiMin), ...
            'obsPhiMax', size(stat.obsPhiMax), ...
            'obsPhiMin', size(stat.obsPhiMin)));
end

function values = make_value_payload(solver)
    values = struct( ...
        'parameters', solver.params, ...
        'p_Psi', solver.p_Psi, ...
        'q_Psi', solver.q_Psi, ...
        'p_var', solver.p_var, ...
        'q_var', solver.q_var);
end

function assembly = make_var_amp_assembly(solver)
    pSize = numel(solver.p_var(:,1));
    qSize = numel(solver.q_var(:,1));
    indexMap = struct();
    indexMap.p_var = reshape(1:pSize,[],1);
    indexMap.q_var = reshape(pSize + (1:qSize),[],1);
    indexMap.varPhiMax = pSize + qSize + 1;
    indexMap.varPhiMin = pSize + qSize + 2;

    assembly = struct( ...
        'indexMap', indexMap, ...
        'parameters', solver.params, ...
        'p_Psi', solver.p_Psi, ...
        'q_Psi', solver.q_Psi, ...
        'p_var', solver.p_var(:,1), ...
        'q_var', solver.q_var(:,1), ...
        'Derivative_obs', [], ...
        'numUnknown', pSize + qSize + 2);
end

function stat = make_guardrail_state(obs, params)
    PV = struct('name', 'var', 'idx', 1);
    t = linspace(0, 2 * pi, 1001).';
    x = [cos(t), sin(t)];
    stat = state(obs, params, t, x, 3, PV);
end

function z = observable_x1(y)
    z = y(:,1);
end

function verifyCallFails(testCase, thunk)
    didFail = false;
    try
        thunk();
    catch
        didFail = true;
    end

    verifyTrue(testCase, didFail);
end
