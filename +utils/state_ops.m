classdef state_ops
    methods (Static)
        function TS_obs = getObs(obs,TS_var)
            TS_obs = fmam_state_ops.getObs(obs,TS_var);
        end

        function [phi_max,phi_min,amplitude,var_max,var_min] = FindExtreme(p,q,L,maxRefinementRounds,extremaResidualTolerance)
            [phi_max,phi_min,amplitude,var_max,var_min] = ...
                fmam_state_ops.FindExtreme( ...
                    p,q,L,maxRefinementRounds,extremaResidualTolerance);
        end

        function output = Trintegration(p,q,phi_L,phi_U)
            output = fmam_state_ops.Trintegration(p,q,phi_L,phi_U);
        end

        function [p,q] = projectFourierSeries(TS,M)
            [p,q] = fmam_state_ops.projectFourierSeries(TS,M);
        end

        function data = buildPrimaryBranchInterpolation(branchX, branchT, branchVar, branchObs)
            data = fmam_state_ops.buildPrimaryBranchInterpolation(branchX, branchT, branchVar, branchObs);
        end

        function [tNorm,TSVarNorm] = normalizePeriodicInputs(t,TS_var)
            [tNorm,TSVarNorm] = fmam_state_ops.normalizePeriodicInputs(t,TS_var);
        end

        function tf = hasRepeatedEndpoint(TS_var)
            tf = fmam_state_ops.hasRepeatedEndpoint(TS_var);
        end

        function issue = timeMapInvariantIssue(PsiValues,p,q,numIntervals)
            issue = fmam_state_ops.timeMapInvariantIssue(PsiValues,p,q,numIntervals);
        end

        function issue = positivePsiIssue(PsiValues)
            issue = fmam_state_ops.positivePsiIssue(PsiValues);
        end

        function assertPositivePsi(PsiValues)
            fmam_state_ops.assertPositivePsi(PsiValues);
        end

        function dt = timeIncrementsFromCoefficients(p,q,numIntervals)
            dt = fmam_state_ops.timeIncrementsFromCoefficients(p,q,numIntervals);
        end

        function issue = positiveTimeIncrementIssue(p,q,numIntervals)
            issue = fmam_state_ops.positiveTimeIncrementIssue(p,q,numIntervals);
        end

        function assertPositiveTimeIncrements(p,q,numIntervals)
            fmam_state_ops.assertPositiveTimeIncrements(p,q,numIntervals);
        end

        function values = evaluateTrigSeries(phi,p,q)
            values = fmam_state_ops.evaluateTrigSeries(phi,p,q);
        end

        function [p_variable,q_variable] = reprojectEqualTimeFourier(tSeries,TS_var,M)
            [p_variable,q_variable] = fmam_state_ops.reprojectEqualTimeFourier(tSeries,TS_var,M);
        end

        function solverView = buildSolverViewFromTrajectory(obs,params,t,TS_var,M,PV,discretization,extremaSearch)
            solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
                obs,params,t,TS_var,M,PV,discretization,extremaSearch);
        end

        function derived = buildDerivedView(obs,solverView,discretization)
            derived = fmam_state_ops.buildDerivedView(obs,solverView,discretization);
        end

        function derived = buildVariableDerivedState(p_var,q_var,p_Psi,q_Psi,Lphi,extremaSearch)
            derived = fmam_state_ops.buildVariableDerivedState( ...
                p_var,q_var,p_Psi,q_Psi,Lphi,extremaSearch);
        end

        function [p_obs,q_obs] = buildObservableFourierCoefficients(obs,p_var,q_var,M,Lphi)
            [p_obs,q_obs] = fmam_state_ops.buildObservableFourierCoefficients(obs,p_var,q_var,M,Lphi);
        end

        function derived = buildObservableDerivedState(obs,p_var,q_var,p_Psi,q_Psi,discretization,extremaSearch)
            derived = fmam_state_ops.buildObservableDerivedState( ...
                obs,p_var,q_var,p_Psi,q_Psi,discretization,extremaSearch);
        end

        function [a,b,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = ...
                reconstructSolverCoefficients(TS_variable,TS_observable,t,M,PV,discretization)
            [a,b,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = ...
                fmam_state_ops.reconstructSolverCoefficients( ...
                    TS_variable,TS_observable,t,M,PV,discretization);
        end

        function value = observableTargetCurrent(obs,p_var,q_var,Phi,k)
            pt_var = fmam_state_ops.evaluateTrigSeries(Phi,p_var,q_var);
            value = obs{k}(pt_var);
        end

        function validatePrimaryVariable(PV,dimVar,dimObs)
            fmam_state_ops.validatePrimaryVariable(PV,dimVar,dimObs);
        end

        function [vc,vs] = vecCS(phi, M, L)
            phi = phi(:);
            if nargin >= 3 && ~isempty(L)
                assert(numel(phi) == L, 'L must match length(phi).');
            else
                L = numel(phi);
            end

            kC = 0:M;
            vc = cos(phi * kC);
            if M > 0
                kS = 1:M;
                vs = sin(phi * kS);
            else
                vs = zeros(L,0,'like',phi);
            end
        end

        function res = residuePhiVar(p,q,phi)
            M = size(q,1);
            [vc,vs] = utils.state_ops.vecCS(phi,M,1);
            res = -vs .* (1:M) * p(2:end,:) + vc(:,2:end) .* (1:M) * q;
        end
    end
end
