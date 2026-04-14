function adapter = state_transaction_adapter(stat)
%STATE_TRANSACTION_ADAPTER Thin transaction boundary around solver-state mutation.
    if ~isa(stat,'state')
        error('state_transaction_adapter:InvalidState', ...
            'stat must be an instance of the class ''state''.')
    end

    adapter = struct();
    adapter.snapshot = @snapshot;
    adapter.restore = @restore;
    adapter.applyUnknownIncrement = @applyUnknownIncrement;
    adapter.currentUnknownInfNorm = @currentUnknownInfNorm;

    function snapshotValue = snapshot()
        snapshotValue = stat.snapshotSolverState();
    end

    function restore(snapshotValue)
        stat.restoreSolverState(snapshotValue);
    end

    function applyUnknownIncrement(meta,increments,scale)
        stat.applyUnknownIncrement(meta,increments,scale);
    end

    function value = currentUnknownInfNorm(unknowns)
        value = stat.currentUnknownInfNorm(unknowns);
    end
end
