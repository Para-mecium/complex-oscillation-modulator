classdef MockFMAMContinuationTask < FMAM_ODE
    properties
        mockFailAboveStep = inf
        mockInvalidPsiAboveStep = inf
        mockIterationsOnSuccess = 3
        mockIterationsOnFailure = 12
    end

    methods
        function obj = MockFMAMContinuationTask(varargin)
            obj@FMAM_ODE(varargin{:});
        end

        function result = fit(obj)
            item = obj.items_perturb(1);
            targetValue = obj.items_per_curr(1);
            currentValue = obj.mockCurrentTargetValue(item);
            stepSize = abs(targetValue - currentValue);

            if stepSize > obj.mockFailAboveStep
                result = obj.mockResult(false,obj.mockIterationsOnFailure,'mock rejection');
                return
            end

            obj.applyMockTarget(item,targetValue);
            if stepSize > obj.mockInvalidPsiAboveStep
                obj.forceInvalidPsi();
            end
            result = obj.mockResult(true,obj.mockIterationsOnSuccess,'mock accepted');
        end
    end

    methods (Access = private)
        function value = mockCurrentTargetValue(obj,item)
            view = obj.exportSolverView();
            switch item.prop
                case 'p_Psi'
                    value = view.p_Psi(item.idx);
                case 'varPhiMax'
                    value = view.varPhiMax(item.idx);
                otherwise
                    error('MockFMAMContinuationTask only supports p_Psi and varPhiMax test targets.')
            end
        end

        function applyMockTarget(obj,item,targetValue)
            view = obj.exportSolverView();
            switch item.prop
                case 'p_Psi'
                    view.p_Psi(item.idx) = targetValue;
                    if item.idx == 1 && view.dimParams >= 1
                        view.params(1) = 1 / targetValue;
                    end
                case 'varPhiMax'
                    view.varPhiMax(item.idx) = targetValue;
                otherwise
                    error('MockFMAMContinuationTask only supports p_Psi and varPhiMax test targets.')
            end
            obj.loadSolverView(view);
        end

        function forceInvalidPsi(obj)
            view = obj.exportSolverView();
            view.p_Psi(:) = 0;
            view.q_Psi(:) = 0;
            obj.loadSolverView(view);
        end

        function result = mockResult(~,converged,iterations,message)
            history = struct( ...
                'iteration', 1, ...
                'objective', 0, ...
                'residualNorm', 0, ...
                'stepNorm', 0, ...
                'acceptedScale', double(converged), ...
                'lambda', 0, ...
                'backtracks', 0, ...
                'accepted', converged, ...
                'maxError', 0, ...
                'solver', 'direct', ...
                'conditionEstimate', 1, ...
                'directConditionEstimate', 1);

            if ~converged
                history.accepted = false;
                history.acceptedScale = 0;
                history.maxError = 1;
            end

            result = struct();
            result.converged = converged;
            result.iterations = iterations;
            result.finalError = zeros(1,4);
            result.message = message;
            result.history = history;
            result.linearResidualNorm = 0;
            result.linearResidual = zeros(1,1);
            result.stepAccepted = converged;
            result.scalarError = double(~converged);
            result.objective = 0;
        end
    end
end
