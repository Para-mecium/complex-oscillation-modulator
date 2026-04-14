classdef ThrowingRestoreState < state
    properties
        throwOnRestoreCall = inf
        restoreCallCount = 0
    end

    methods
        function obj = ThrowingRestoreState(varargin)
            obj@state(varargin{:});
        end

        function restoreSolverState(obj,snapshot)
            obj.restoreCallCount = obj.restoreCallCount + 1;
            if obj.restoreCallCount == obj.throwOnRestoreCall
                error('ThrowingRestoreState:RestoreFailure', ...
                    'Injected restore failure during trial activation.');
            end

            restoreSolverState@state(obj,snapshot);
        end
    end
end
