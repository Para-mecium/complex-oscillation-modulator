function [A, res, indexMap, unknowns] = assemble_newton_linear_system(ctx)
%ASSEMBLE_NEWTON_LINEAR_SYSTEM Build the Newton matrix and residual vector.

    p_Psi_init_ = ctx.p_Psi_init;
    q_Psi_init_ = ctx.q_Psi_init;

    system = ctx.sys;
    observable = ctx.obs;
    Derivative_var = ctx.derivatives.var;
    Derivative_obs = ctx.derivatives.obs;
    DDerivative_obs = ctx.derivatives.obs2;
    L = ctx.discretization.assemblySampleCount;
    N = ctx.dimVar;
    n = ctx.dimObs;
    MVar = ctx.truncationOrder;
    m = ctx.dimParams;

    parameters = ctx.params;
    Coe_Controlled = ctx.items_controlled;
    if isfield(ctx, 'targetRuleCtx') && ~isempty(ctx.targetRuleCtx)
        targetCtx = ctx.targetRuleCtx;
    else
        targetCtx = build_target_rules_context(ctx);
    end

    p_var = ctx.p_var;
    q_var = ctx.q_var;
    p_Psi = ctx.p_Psi;
    q_Psi = ctx.q_Psi;
    MPsi = size(q_Psi,1);

    phi = (0:L-1)' * 2 * pi / L;
    [vcVar, vsVar] = FMAM_ODE.Vec_CS(phi, MVar, L);
    [vcPsi, vsPsi] = FMAM_ODE.Vec_CS(phi, MPsi, L);
    TS_var = vcVar * p_var + vsVar * q_var;
    Psi = vcPsi * p_Psi + vsPsi * q_Psi;
    TS_Dvar = -(vsVar .* (1:MVar)) * p_var(2:end, :) + (vcVar(:, 2:end) .* (1:MVar)) * q_var;

    unknowns = {'params', 'p_Psi', 'q_Psi', 'p_var', 'q_var', ...
        'varPhiMax', 'varPhiMin', 'obsPhiMax', 'obsPhiMin'};

    indexMap = struct();
    lastIdx = 0;
    for i = 1:numel(unknowns)
        propname = unknowns{i};
        rows = size(ctx.(propname), 1);
        columns = size(ctx.(propname), 2);
        indexMap.(propname) = reshape(lastIdx + (1:rows * columns), [], columns);
        lastIdx = lastIdx + rows * columns;
    end
    A = zeros(lastIdx, lastIdx);
    res = zeros(lastIdx, 1);

    idxCollum = 1:max([indexMap.params(:); indexMap.p_Psi(:); indexMap.q_Psi(:); ...
        indexMap.p_var(:); indexMap.q_var(:)]);
    idxLastRow = 0;
    for i = 1:N
        [coe_params, coe_p_Psi, coe_q_Psi, coe_p_var, coe_q_var] = ...
            FMAM_ODE.delta_coe_system(system, Derivative_var, Psi, TS_var, ...
            parameters, vcPsi, vsPsi, vcVar, vsVar, i);
        A_temp = fft([coe_params, coe_p_Psi, coe_q_Psi, coe_p_var, coe_q_var]) / L;
        A_plus1 = real(A_temp);
        A_plus2 = imag(A_temp);

        residue = FMAM_ODE.residue_system(system, Psi, TS_var, TS_Dvar, parameters, i);
        res_temp = fft(residue) / L;
        res_plus1 = real(res_temp);
        res_plus2 = imag(res_temp);

        A(idxLastRow + (1:MVar + 1), idxCollum) = A_plus1(1:MVar + 1, :);
        res(idxLastRow + (1:MVar + 1)) = res_plus1(1:MVar + 1);
        A(idxLastRow + MVar + 1 + (1:MVar), idxCollum) = A_plus2(2:MVar + 1, :);
        res(idxLastRow + MVar + 1 + (1:MVar)) = res_plus2(2:MVar + 1);

        idxLastRow = idxLastRow + (2 * MVar + 1);

        if ~ctx.isPsiUpdated && should_insert_frozen_psi_closure_row(ctx.PV, i)
            A(idxLastRow + 1, idxCollum) = A_plus1(MVar + 2, :);
            res(idxLastRow + 1) = res_plus1(MVar + 2);
            idxLastRow = idxLastRow + 1;
        end
    end

    for i = 1:N
        [needMax, needMin] = fmam_target_rules('needs_variable_extrema', targetCtx, ctx.items_perturb, i);
        idx_phi_max = indexMap.varPhiMax(i);
        idx_phi_min = indexMap.varPhiMin(i);
        if needMax
            p = p_var(:, i);
            q = q_var(:, i);
            idx_p_var_i = indexMap.p_var(:, i);
            idx_q_var_i = indexMap.q_var(:, i);
            phi = ctx.varPhiMax(i);

            [coe_p_var_i, coe_q_var_i, coe_phi] = FMAM_ODE.delta_coe_phi_var(p, q, phi);

            A(idxLastRow + 1, idx_p_var_i) = coe_p_var_i;
            A(idxLastRow + 1, idx_q_var_i) = coe_q_var_i;
            A(idxLastRow + 1, idx_phi_max) = coe_phi;
            res(idxLastRow + 1) = -FMAM_ODE.residue_phi_var(p, q, phi);
        else
            A(idxLastRow + 1, idx_phi_max) = 1;
            res(idxLastRow + 1) = 0;
        end
        idxLastRow = idxLastRow + 1;

        if needMin
            p = p_var(:, i);
            q = q_var(:, i);
            idx_p_var_i = indexMap.p_var(:, i);
            idx_q_var_i = indexMap.q_var(:, i);
            phi = ctx.varPhiMin(i);

            [coe_p_var_i, coe_q_var_i, coe_phi] = FMAM_ODE.delta_coe_phi_var(p, q, phi);

            A(idxLastRow + 1, idx_p_var_i) = coe_p_var_i;
            A(idxLastRow + 1, idx_q_var_i) = coe_q_var_i;
            A(idxLastRow + 1, idx_phi_min) = coe_phi;
            res(idxLastRow + 1) = -FMAM_ODE.residue_phi_var(p, q, phi);
        else
            A(idxLastRow + 1, idx_phi_min) = 1;
            res(idxLastRow + 1) = 0;
        end
        idxLastRow = idxLastRow + 1;
    end

    for k = 1:n
        [needMax, needMin] = fmam_target_rules('needs_observable_extrema', targetCtx, ctx.items_perturb, k);
        idx_phi_max = indexMap.obsPhiMax(k);
        idx_phi_min = indexMap.obsPhiMin(k);

        if needMax
            phi = ctx.obsPhiMax(k);
            [coe_p_var, coe_q_var, coe_phi] = ...
                FMAM_ODE.delta_coe_obs_phi(Derivative_obs, DDerivative_obs, p_var, q_var, phi, k);

            for i = 1:N
                idx_p_var_i = indexMap.p_var(:, i);
                idx_q_var_i = indexMap.q_var(:, i);
                A(idxLastRow + 1, idx_p_var_i) = coe_p_var(:, i)';
                A(idxLastRow + 1, idx_q_var_i) = coe_q_var(:, i)';
            end
            A(idxLastRow + 1, idx_phi_max) = coe_phi;
            res(idxLastRow + 1) = -FMAM_ODE.residue_phi_obs(Derivative_obs, p_var, q_var, phi, k);
        else
            A(idxLastRow + 1, idx_phi_max) = 1;
            res(idxLastRow + 1) = 0;
        end
        idxLastRow = idxLastRow + 1;

        if needMin
            phi = ctx.obsPhiMin(k);
            [coe_p_var, coe_q_var, coe_phi] = ...
                FMAM_ODE.delta_coe_obs_phi(Derivative_obs, DDerivative_obs, p_var, q_var, phi, k);

            for i = 1:N
                idx_p_var_i = indexMap.p_var(:, i);
                idx_q_var_i = indexMap.q_var(:, i);
                A(idxLastRow + 1, idx_p_var_i) = coe_p_var(:, i)';
                A(idxLastRow + 1, idx_q_var_i) = coe_q_var(:, i)';
            end
            A(idxLastRow + 1, idx_phi_min) = coe_phi;
            res(idxLastRow + 1) = -FMAM_ODE.residue_phi_obs(Derivative_obs, p_var, q_var, phi, k);
        else
            A(idxLastRow + 1, idx_phi_min) = 1;
            res(idxLastRow + 1) = 0;
        end
        idxLastRow = idxLastRow + 1;
    end

    idx_params_fix = 1:m;
    idx_params_fix(Coe_Controlled) = [];
    for i = idx_params_fix
        A(idxLastRow + 1, i) = 1;
        res(idxLastRow + 1) = 0;
        idxLastRow = idxLastRow + 1;
    end

    assembly = build_target_row_assembly(indexMap, parameters, p_Psi, q_Psi, ...
        p_var, q_var, Derivative_obs, size(A, 2));
    for j = 1:numel(ctx.items_perturb)
        item_per = ctx.items_perturb(j);
        [targetRow, targetResidual] = fmam_target_rules('target_row', ...
            targetCtx, item_per, ctx.targetCurr(j), assembly);
        A(idxLastRow + 1, :) = targetRow;
        res(idxLastRow + 1) = targetResidual;
        idxLastRow = idxLastRow + 1;
    end

    if ctx.isPsiUpdated
        PV = ctx.PV;
        if strcmpi(PV.name, 'var')
            idx_PV = PV.idx;
            idx_p = indexMap.p_var(3:end, idx_PV);
            idx_q = indexMap.q_var(:, idx_PV);

            if ~isempty(idx_p)
                A(idxLastRow + (1:MVar - 1), idx_p) = diag(ones(MVar - 1, 1));
                res(idxLastRow + (1:MVar - 1)) = -p_var(3:end, idx_PV);
                idxLastRow = idxLastRow + MVar - 1;
            end

            if ~isempty(idx_q)
                A(idxLastRow + (1:MVar), idx_q) = diag(ones(MVar, 1));
                res(idxLastRow + (1:MVar)) = -q_var(:, idx_PV);
                idxLastRow = idxLastRow + MVar;
            end
        elseif strcmpi(PV.name, 'obs')
            [coe_p_var, coe_q_var] = FMAM_ODE.delta_coe_observable(Derivative_obs, p_var, q_var, PV, L);
            idx_column = [indexMap.p_var(:); indexMap.q_var(:)]';
            A_temp = fft([coe_p_var, coe_q_var]) / L;
            A_plus1 = real(A_temp);
            A_plus2 = imag(A_temp);

            f = observable{PV.idx};
            res_obs = f(TS_var);
            res_temp = fft(res_obs) / L;
            res_plus1 = real(res_temp);
            res_plus2 = imag(res_temp);

            if MVar > 1
                A(idxLastRow + (1:MVar - 1), idx_column) = A_plus1(3:MVar + 1, :);
                res(idxLastRow + (1:MVar - 1)) = -res_plus1(3:MVar + 1);
                idxLastRow = idxLastRow + MVar - 1;
            end

            A(idxLastRow + (1:MVar), idx_column) = A_plus2(2:MVar + 1, :);
            res(idxLastRow + (1:MVar)) = -res_plus2(2:MVar + 1);
            idxLastRow = idxLastRow + MVar;
        end
    else
        idx_p_Psi = indexMap.p_Psi(2:end);
        idx_q_Psi = indexMap.q_Psi;
        if ~isempty(idx_p_Psi)
            numFrozenCos = numel(idx_p_Psi);
            A(idxLastRow + (1:numFrozenCos), idx_p_Psi) = diag(ones(numFrozenCos, 1));
            res(idxLastRow + (1:numFrozenCos)) = -(p_Psi(2:end) - p_Psi_init_(2:end));
            idxLastRow = idxLastRow + numFrozenCos;
        end

        if ~isempty(idx_q_Psi)
            numFrozenSin = numel(idx_q_Psi);
            A(idxLastRow + (1:numFrozenSin), idx_q_Psi) = diag(ones(numFrozenSin, 1));
            res(idxLastRow + (1:numFrozenSin)) = -(q_Psi - q_Psi_init_);
            idxLastRow = idxLastRow + numFrozenSin;
        end
    end

    if idxLastRow ~= size(A, 1)
        error('FMAM_ODE:AssemblySizeMismatch', ...
            'Newton assembly produced %d rows for %d unknowns.', idxLastRow, size(A, 1));
    end
end

function tf = should_insert_frozen_psi_closure_row(PV, eqIdx)
    if strcmpi(PV.name, 'var')
        tf = (eqIdx == PV.idx);
    else
        tf = (eqIdx == 1);
    end
end

function targetCtx = build_target_rules_context(ctx)
    discretization = fmam_state_defaults.normalizeDiscretization(ctx.discretization);
    extremaSearch = fmam_state_ops.normalizeExtremaSearchSettings(ctx.extremaSearch);
    targetCtx = struct();
    targetCtx.obs = ctx.obs;
    targetCtx.dimVar = ctx.dimVar;
    targetCtx.dimObs = ctx.dimObs;
    targetCtx.dimParams = ctx.dimParams;
    targetCtx.propertySizes = struct( ...
        'p_Psi', size(ctx.p_Psi), ...
        'q_Psi', size(ctx.q_Psi), ...
        'p_var', size(ctx.p_var), ...
        'q_var', size(ctx.q_var), ...
        'varPhiMax', size(ctx.varPhiMax), ...
        'varPhiMin', size(ctx.varPhiMin), ...
        'obsPhiMax', size(ctx.obsPhiMax), ...
        'obsPhiMin', size(ctx.obsPhiMin), ...
        'params', size(ctx.params));
    targetCtx.solver = struct( ...
        'params', ctx.params, ...
        'p_Psi', ctx.p_Psi, ...
        'q_Psi', ctx.q_Psi, ...
        'p_var', ctx.p_var, ...
        'q_var', ctx.q_var, ...
        'varPhiMax', ctx.varPhiMax, ...
        'varPhiMin', ctx.varPhiMin, ...
        'obsPhiMax', ctx.obsPhiMax, ...
        'obsPhiMin', ctx.obsPhiMin, ...
        'PV', ctx.PV, ...
        'propertySizes', targetCtx.propertySizes);
    targetCtx.discretization = discretization;
    targetCtx.extremaSearch = extremaSearch;
    targetCtx.derived = fmam_state_ops.buildDerivedView( ...
        ctx.obs,targetCtx.solver,discretization);
end

function assembly = build_target_row_assembly(indexMap, parameters, p_Psi, q_Psi, ...
        p_var, q_var, Derivative_obs, numUnknown)
    assembly = struct();
    assembly.indexMap = indexMap;
    assembly.parameters = parameters;
    assembly.p_Psi = p_Psi;
    assembly.q_Psi = q_Psi;
    assembly.p_var = p_var;
    assembly.q_var = q_var;
    assembly.Derivative_obs = Derivative_obs;
    assembly.numUnknown = numUnknown;
end
