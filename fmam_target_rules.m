function varargout = fmam_target_rules(action,ctx,varargin)
%FMAM_TARGET_RULES Shared target semantics for FMAM modulation targets.
% Supported actions:
%   'validate_item'           -> item = fmam_target_rules(action, ctx, item)
%   'needs_variable_extrema'  -> [needMax,needMin] = fmam_target_rules(action, ctx, items, varIdx)
%   'needs_observable_extrema'-> [needMax,needMin] = fmam_target_rules(action, ctx, items, obsIdx)
%   'current_value'           -> value = fmam_target_rules(action, ctx, item, values)
%   'linear_index'            -> idx = fmam_target_rules(action, ctx, propname, rawIdx)
%   'continuation_delta'      -> [delta,isWrapped,period] = fmam_target_rules(action, ctx, item, startValue, rawTarget)
%   'target_row'              -> [row,residual] = fmam_target_rules(action, ctx, item, targetValue, assembly)

    switch lower(action)
        case 'validate_item'
            varargout{1} = validate_item(ctx,varargin{1});
        case 'needs_variable_extrema'
            [varargout{1},varargout{2}] = needs_variable_extrema(varargin{1},varargin{2});
        case 'needs_observable_extrema'
            [varargout{1},varargout{2}] = needs_observable_extrema(varargin{1},varargin{2});
        case 'current_value'
            varargout{1} = current_value(ctx,varargin{1},varargin{2});
        case 'linear_index'
            varargout{1} = linear_index(ctx,varargin{1},varargin{2});
        case 'continuation_delta'
            [varargout{1},varargout{2},varargout{3}] = continuation_delta(ctx,varargin{1},varargin{2},varargin{3});
        case 'target_row'
            [varargout{1},varargout{2}] = target_row(ctx,varargin{1},varargin{2},varargin{3});
        otherwise
            error('fmam_target_rules:UnknownAction', ...
                'Unsupported target rules action ''%s''.',action)
    end
end

function item = validate_item(ctx,item)
    item.prop = fmam_target_canonicalize(item.prop);
    idx = item.idx;
    if ~isnumeric(idx) || isempty(idx) || any(~isfinite(idx)) || any(idx ~= floor(idx))
        error('Target indices for %s must be finite integers.',item.prop)
    end

    switch item.prop
        case {'params','p_Psi','q_Psi','varPhiMax','varPhiMin','obsPhiMax','obsPhiMin', ...
                'varAmp','varMax','varMin','obsAmp','obsMax','obsMin'}
            if numel(idx) ~= 1
                error('%s targets require a single scalar index.',item.prop)
            end
        case {'p_var','q_var','varPhase'}
            if numel(idx) ~= 2
                error('%s targets require a two-component index.',item.prop)
            end
        otherwise
            error('Unsupported target property ''%s''.',item.prop)
    end

    switch item.prop
        case 'params'
            upperBound = ctx.dimParams;
        case 'p_Psi'
            upperBound = ctx.propertySizes.p_Psi(1);
        case 'q_Psi'
            upperBound = ctx.propertySizes.q_Psi(1);
        case {'varPhiMax','varPhiMin','varAmp','varMax','varMin'}
            upperBound = ctx.dimVar;
        case {'obsPhiMax','obsPhiMin','obsAmp','obsMax','obsMin'}
            upperBound = ctx.dimObs;
            if upperBound == 0
                error('%s targets require at least one observable.',item.prop)
            end
        case 'p_var'
            validate_matrix_index(idx,ctx.propertySizes.p_var,'p_var index exceeds the size of stat.p_var.')
            return
        case 'q_var'
            validate_matrix_index(idx,ctx.propertySizes.q_var,'q_var index exceeds the size of stat.q_var.')
            return
        case 'varPhase'
            upperBound = ctx.dimVar;
        otherwise
            error('Unsupported target property ''%s''.',item.prop)
    end

    if any(idx < 1) || any(idx > upperBound)
        error('Index for %s exceeds the supported range.',item.prop)
    end
end

function [needMax,needMin] = needs_variable_extrema(items,varIdx)
    needMax = false;
    needMin = false;
    for j = 1:numel(items)
        item = items(j);
        switch item.prop
            case 'varAmp'
                if item.idx == varIdx
                    needMax = true;
                    needMin = true;
                end
            case {'varPhiMax','varMax'}
                if item.idx == varIdx
                    needMax = true;
                end
            case {'varPhiMin','varMin'}
                if item.idx == varIdx
                    needMin = true;
                end
            case 'varPhase'
                if any(item.idx == varIdx)
                    needMax = true;
                end
        end

        if needMax && needMin
            break
        end
    end
end

function [needMax,needMin] = needs_observable_extrema(items,obsIdx)
    needMax = false;
    needMin = false;
    for j = 1:numel(items)
        item = items(j);
        switch item.prop
            case 'obsAmp'
                if item.idx == obsIdx
                    needMax = true;
                    needMin = true;
                end
            case {'obsPhiMax','obsMax'}
                if item.idx == obsIdx
                    needMax = true;
                end
            case {'obsPhiMin','obsMin'}
                if item.idx == obsIdx
                    needMin = true;
                end
        end

        if needMax && needMin
            break
        end
    end
end

function value = current_value(ctx,item,values)
    switch item.prop
        case 'params'
            if isfield(values,'parameters')
                value = values.parameters(item.idx);
            else
                value = values.params(item.idx);
            end
        case 'p_Psi'
            value = values.p_Psi(item.idx);
        case 'q_Psi'
            value = values.q_Psi(item.idx);
        case 'p_var'
            value = values.p_var(item.idx(1),item.idx(2));
        case 'q_var'
            value = values.q_var(item.idx(1),item.idx(2));
        case 'varPhiMax'
            value = ctx.solver.varPhiMax(item.idx);
        case 'varPhiMin'
            value = ctx.solver.varPhiMin(item.idx);
        case 'obsPhiMax'
            value = ctx.solver.obsPhiMax(item.idx);
        case 'obsPhiMin'
            value = ctx.solver.obsPhiMin(item.idx);
        case 'varAmp'
            value = ctx.derived.varAmp(item.idx);
        case 'varMax'
            value = ctx.derived.varMax(item.idx);
        case 'varMin'
            value = ctx.derived.varMin(item.idx);
        case 'obsAmp'
            value = ctx.derived.obsAmp(item.idx);
        case 'obsMax'
            value = ctx.derived.obsMax(item.idx);
        case 'obsMin'
            value = ctx.derived.obsMin(item.idx);
        case 'varPhase'
            value = ctx.derived.varPhase(item.idx(1),item.idx(2));
        otherwise
            error('Unsupported modulation target property ''%s''.',item.prop)
    end
end

function idx = linear_index(ctx,propname,rawIdx)
    if isscalar(rawIdx)
        idx = rawIdx;
        return
    end

    idxCell = num2cell(rawIdx);
    idx = sub2ind(ctx.propertySizes.(propname),idxCell{:});
end

function [deltaValue,isWrapped,wrapPeriod] = continuation_delta(ctx,item,startValue,rawTarget)
    isWrapped = false;
    wrapPeriod = NaN;
    deltaValue = rawTarget - startValue;

    switch item.prop
        case {'varPhiMax','varPhiMin','obsPhiMax','obsPhiMin'}
            isWrapped = true;
            wrapPeriod = 2 * pi;
            deltaValue = wrap_periodic_difference(deltaValue,wrapPeriod);
        case 'varPhase'
            isWrapped = true;
            wrapPeriod = max(abs(ctx.derived.period),eps);
            deltaValue = wrap_periodic_difference(deltaValue,wrapPeriod);
    end
end

function [row,residual] = target_row(ctx,item,targetValue,assembly)
    row = zeros(1,assembly.numUnknown);
    values = struct( ...
        'parameters', assembly.parameters, ...
        'p_Psi', assembly.p_Psi, ...
        'q_Psi', assembly.q_Psi, ...
        'p_var', assembly.p_var, ...
        'q_var', assembly.q_var);
    currentValue = current_value(ctx,item,values);

    switch item.prop
        case {'params','p_Psi','q_Psi','p_var','q_var','varPhiMax','varPhiMin','obsPhiMax','obsPhiMin'}
            idxTarget = linear_index(ctx,item.prop,item.idx);
            row(assembly.indexMap.(item.prop)(idxTarget)) = 1;

        case 'varAmp'
            idxTarget = item.idx;
            p = assembly.p_var(:,idxTarget);
            q = assembly.q_var(:,idxTarget);
            phiMax = ctx.solver.varPhiMax(idxTarget);
            phiMin = ctx.solver.varPhiMin(idxTarget);

            [coe_p_var_max,coe_q_var_max,coe_phi_var_max] = FMAM_ODE.delta_coe_var_target(p,q,phiMax);
            [coe_p_var_min,coe_q_var_min,coe_phi_var_min] = FMAM_ODE.delta_coe_var_target(p,q,phiMin);

            row(assembly.indexMap.p_var(:,idxTarget)) = 0.5 * (coe_p_var_max - coe_p_var_min);
            row(assembly.indexMap.q_var(:,idxTarget)) = 0.5 * (coe_q_var_max - coe_q_var_min);
            row(assembly.indexMap.varPhiMax(idxTarget)) = 0.5 * coe_phi_var_max;
            row(assembly.indexMap.varPhiMin(idxTarget)) = -0.5 * coe_phi_var_min;

        case 'varMax'
            idxTarget = item.idx;
            [coe_p_var,coe_q_var,coe_phi_var] = FMAM_ODE.delta_coe_var_target(assembly.p_var(:,idxTarget),assembly.q_var(:,idxTarget),ctx.solver.varPhiMax(idxTarget));
            row(assembly.indexMap.p_var(:,idxTarget)) = coe_p_var;
            row(assembly.indexMap.q_var(:,idxTarget)) = coe_q_var;
            row(assembly.indexMap.varPhiMax(idxTarget)) = coe_phi_var;

        case 'varMin'
            idxTarget = item.idx;
            [coe_p_var,coe_q_var,coe_phi_var] = FMAM_ODE.delta_coe_var_target(assembly.p_var(:,idxTarget),assembly.q_var(:,idxTarget),ctx.solver.varPhiMin(idxTarget));
            row(assembly.indexMap.p_var(:,idxTarget)) = coe_p_var;
            row(assembly.indexMap.q_var(:,idxTarget)) = coe_q_var;
            row(assembly.indexMap.varPhiMin(idxTarget)) = coe_phi_var;

        case 'obsAmp'
            idxTarget = item.idx;
            phiMax = ctx.solver.obsPhiMax(idxTarget);
            phiMin = ctx.solver.obsPhiMin(idxTarget);
            [coe_p_var_max,coe_q_var_max,coe_phi_obs_max] = FMAM_ODE.delta_coe_obs_target(assembly.Derivative_obs,assembly.p_var,assembly.q_var,phiMax,idxTarget);
            [coe_p_var_min,coe_q_var_min,coe_phi_obs_min] = FMAM_ODE.delta_coe_obs_target(assembly.Derivative_obs,assembly.p_var,assembly.q_var,phiMin,idxTarget);

            for i = 1:ctx.dimVar
                row(assembly.indexMap.p_var(:,i)) = 0.5 * (coe_p_var_max(:,i)' - coe_p_var_min(:,i)');
                row(assembly.indexMap.q_var(:,i)) = 0.5 * (coe_q_var_max(:,i)' - coe_q_var_min(:,i)');
            end
            row(assembly.indexMap.obsPhiMax(idxTarget)) = 0.5 * coe_phi_obs_max;
            row(assembly.indexMap.obsPhiMin(idxTarget)) = -0.5 * coe_phi_obs_min;

        case 'obsMax'
            idxTarget = item.idx;
            [coe_p_var,coe_q_var,coe_phi_obs] = FMAM_ODE.delta_coe_obs_target(assembly.Derivative_obs,assembly.p_var,assembly.q_var,ctx.solver.obsPhiMax(idxTarget),idxTarget);
            for i = 1:ctx.dimVar
                row(assembly.indexMap.p_var(:,i)) = coe_p_var(:,i)';
                row(assembly.indexMap.q_var(:,i)) = coe_q_var(:,i)';
            end
            row(assembly.indexMap.obsPhiMax(idxTarget)) = coe_phi_obs;

        case 'obsMin'
            idxTarget = item.idx;
            [coe_p_var,coe_q_var,coe_phi_obs] = FMAM_ODE.delta_coe_obs_target(assembly.Derivative_obs,assembly.p_var,assembly.q_var,ctx.solver.obsPhiMin(idxTarget),idxTarget);
            for i = 1:ctx.dimVar
                row(assembly.indexMap.p_var(:,i)) = coe_p_var(:,i)';
                row(assembly.indexMap.q_var(:,i)) = coe_q_var(:,i)';
            end
            row(assembly.indexMap.obsPhiMin(idxTarget)) = coe_phi_obs;

        case 'varPhase'
            idxVar1 = item.idx(1);
            idxVar2 = item.idx(2);
            [coe_p_Psi,coe_q_Psi,coe_phi1,coe_phi2] = FMAM_ODE.delta_coe_state_phase(assembly.p_Psi,assembly.q_Psi,ctx.solver.varPhiMax(idxVar1),ctx.solver.varPhiMax(idxVar2));

            row(assembly.indexMap.p_Psi) = coe_p_Psi;
            row(assembly.indexMap.q_Psi) = coe_q_Psi;
            row(assembly.indexMap.varPhiMax(idxVar1)) = coe_phi1;
            row(assembly.indexMap.varPhiMax(idxVar2)) = coe_phi2;

        otherwise
            error('Unsupported modulation target property ''%s''.',item.prop)
    end

    residual = targetValue - currentValue;
end

function validate_matrix_index(idx,bounds,messageText)
    if any(idx(:)' < 1) || idx(1) > bounds(1) || idx(2) > bounds(2)
        error(messageText)
    end
end

function delta = wrap_periodic_difference(delta,period)
    if ~(isfinite(period) && period > 0)
        return
    end
    halfPeriod = period / 2;
    delta = mod(delta + halfPeriod,period) - halfPeriod;
end
