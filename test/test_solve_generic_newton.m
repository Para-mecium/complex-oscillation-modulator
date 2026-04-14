function tests = test_solve_generic_newton
%TEST_SOLVE_GENERIC_NEWTON Regression coverage for the generic Newton solver.

    tests = functiontests(localfunctions);
end

function testWellConditionedSystemPrefersDirectSolve(testCase)
    x = [0; 0];
    J = [2, 0; 0, 1];
    target = [2; 1];

    problem = struct();
    problem.linearize = @linearize;
    problem.snapshot = @snapshot;
    problem.restore = @restore;
    problem.applyIncrement = @applyIncrement;
    problem.measure = @measure;

    opts = struct('maxIterations', 4, 'incrementTolerance', 1e-12, ...
        'initialLambda', 1e-6, 'lambdaMin', 1e-12, 'lambdaMax', 1e6, ...
        'lambdaGrow', 10, 'lambdaShrink', 0.3, ...
        'directConditionThreshold', 1e-10, 'lmConditionThreshold', 1e-12);

    result = solve_generic_newton(problem, opts);

    verifyTrue(testCase, result.converged);
    verifyTrue(testCase, result.stepAccepted);
    verifyGreaterThanOrEqual(testCase, numel(result.history), 1);
    verifyTrue(testCase, all([result.history.accepted]));
    verifyTrue(testCase, all(strcmp({result.history.solver}, 'direct')));
    verifyEqual(testCase, [result.history.lambda], zeros(1, numel(result.history)), 'AbsTol', 0);
    verifyGreaterThan(testCase, min([result.history.conditionEstimate]), opts.directConditionThreshold);
    verifyLessThan(testCase, norm(target - J * x, inf), 1e-10);

    function [JCurrent, residual, meta] = linearize()
        JCurrent = J;
        residual = target - J * x;
        meta = struct('iterateNorm', norm(x, inf));
    end

    function value = snapshot()
        value = x;
    end

    function restore(snapshot)
        x = snapshot;
    end

    function applyIncrement(~, delta, scale)
        x = x + scale * delta;
    end

    function measurement = measure()
        err = norm(target - J * x, inf);
        measurement = struct('errorVec', err, 'converged', err <= 1e-10, 'scalarError', err);
    end
end

function testIllConditionedSystemUsesLM(testCase)
    x = [0; 0];
    J = [1, 0; 1, 1e-12; 0, 1e-12];
    target = [1; 1; 0];

    problem = struct();
    problem.linearize = @linearize;
    problem.snapshot = @snapshot;
    problem.restore = @restore;
    problem.applyIncrement = @applyIncrement;
    problem.measure = @measure;

    opts = struct('maxIterations', 6, 'incrementTolerance', 1e-12, ...
        'initialLambda', 1e-8, 'lambdaMin', 1e-12, 'lambdaMax', 1e6, ...
        'lambdaGrow', 10, 'lambdaShrink', 0.3, ...
        'directConditionThreshold', 1e-10, 'lmConditionThreshold', 1e-12);

    result = solve_generic_newton(problem, opts);

    verifyTrue(testCase, result.converged);
    verifyTrue(testCase, result.stepAccepted);
    verifyGreaterThanOrEqual(testCase, numel(result.history), 1);
    verifyTrue(testCase, all(strcmp({result.history.solver}, 'lm')));
    verifyTrue(testCase, all([result.history.lambda] > 0));
    verifyLessThan(testCase, norm(target - J * x, inf), 1e-8);

    function [JCurrent, residual, meta] = linearize()
        JCurrent = J;
        residual = target - J * x;
        meta = struct('iterateNorm', norm(x, inf));
    end

    function value = snapshot()
        value = x;
    end

    function restore(snapshot)
        x = snapshot;
    end

    function applyIncrement(~, delta, scale)
        x = x + scale * delta;
    end

    function measurement = measure()
        err = norm(target - J * x, inf);
        measurement = struct('errorVec', err, 'converged', err <= 1e-8, 'scalarError', err);
    end
end

function testFailureRestoresStateAfterRejectedCandidates(testCase)
    x = 0;

    problem = struct();
    problem.linearize = @linearize;
    problem.snapshot = @snapshot;
    problem.restore = @restore;
    problem.applyIncrement = @applyIncrement;
    problem.measure = @measure;

    opts = struct('maxIterations', 1, 'initialLambda', 1e-8, 'lambdaMin', 1e-8, ...
        'lambdaMax', 1e-4, 'lambdaGrow', 10, 'lambdaShrink', 0.3, ...
        'directConditionThreshold', 1e-10, 'lmConditionThreshold', 1e-12);

    result = solve_generic_newton(problem, opts);

    verifyFalse(testCase, result.converged);
    verifyFalse(testCase, result.stepAccepted);
    verifyEqual(testCase, x, 0, 'AbsTol', 0);
    verifyEqual(testCase, result.errorVec, 1, 'AbsTol', 0);

    function [JCurrent, residual, meta] = linearize()
        JCurrent = 1;
        residual = 1 - x;
        meta = struct('iterateNorm', abs(x));
    end

    function value = snapshot()
        value = x;
    end

    function restore(snapshot)
        x = snapshot;
    end

    function applyIncrement(~, delta, scale)
        x = x + scale * delta;
        error('toy:RejectedCandidate', 'Reject every candidate after mutating state.');
    end

    function measurement = measure()
        err = abs(1 - x);
        measurement = struct('errorVec', err, 'converged', false, 'scalarError', err);
    end
end

function testCandidateValidatorBacktracksAndAcceptsScaledStep(testCase)
    x = 0;

    problem = struct();
    problem.linearize = @linearize;
    problem.snapshot = @snapshot;
    problem.restore = @restore;
    problem.applyIncrement = @applyIncrement;
    problem.measure = @measure;
    problem.validateCandidate = @validateCandidate;

    opts = struct('maxIterations', 2, 'incrementTolerance', 1e-12, ...
        'initialLambda', 1e-8, 'lambdaMin', 1e-12, 'lambdaMax', 1e6, ...
        'lambdaGrow', 10, 'lambdaShrink', 0.3, ...
        'directConditionThreshold', 1e-10, 'lmConditionThreshold', 1e-12, ...
        'candidateBacktrackingFactor', 0.5, ...
        'candidateBacktrackingMaxBacktracks', 4);

    result = solve_generic_newton(problem, opts);

    verifyTrue(testCase, result.converged);
    verifyTrue(testCase, result.stepAccepted);
    verifyEqual(testCase, numel(result.history), 1);
    verifyEqual(testCase, result.history.acceptedScale, 0.5, 'AbsTol', 0);
    verifyEqual(testCase, result.history.backtracks, 1);
    verifyEqual(testCase, x, 0.5, 'AbsTol', 1e-12);

    function [JCurrent, residual, meta] = linearize()
        JCurrent = 1;
        residual = 1 - x;
        meta = struct('iterateNorm', abs(x));
    end

    function value = snapshot()
        value = x;
    end

    function restore(snapshot)
        x = snapshot;
    end

    function applyIncrement(~, delta, scale)
        x = x + scale * delta;
    end

    function measurement = measure()
        err = abs(0.5 - x);
        measurement = struct('errorVec', err, 'converged', err <= 1e-12, 'scalarError', err);
    end

    function validation = validateCandidate(~, ~, ~)
        validation = struct('isValid', x <= 0.5, 'message', 'x must stay below 0.5');
    end
end

function testCandidateValidatorStopsAtScaledIncrementTolerance(testCase)
    x = 0;

    problem = struct();
    problem.linearize = @linearize;
    problem.snapshot = @snapshot;
    problem.restore = @restore;
    problem.applyIncrement = @applyIncrement;
    problem.measure = @measure;
    problem.validateCandidate = @validateCandidate;

    opts = struct('maxIterations', 1, 'incrementTolerance', 0.2, ...
        'initialLambda', 1e-8, 'lambdaMin', 1e-12, 'lambdaMax', 1e6, ...
        'lambdaGrow', 10, 'lambdaShrink', 0.3, ...
        'directConditionThreshold', 1e-10, 'lmConditionThreshold', 1e-12, ...
        'candidateBacktrackingFactor', 0.5, ...
        'candidateBacktrackingMaxBacktracks', 8);

    result = solve_generic_newton(problem, opts);

    verifyFalse(testCase, result.converged);
    verifyFalse(testCase, result.stepAccepted);
    verifyEqual(testCase, result.message, 'increment tolerance satisfied');
    verifyEqual(testCase, x, 0, 'AbsTol', 0);
    verifyEqual(testCase, result.history.backtracks, 3);
    verifyEqual(testCase, result.history.acceptedScale, 0, 'AbsTol', 0);

    function [JCurrent, residual, meta] = linearize()
        JCurrent = 1;
        residual = 1 - x;
        meta = struct('iterateNorm', abs(x));
    end

    function value = snapshot()
        value = x;
    end

    function restore(snapshot)
        x = snapshot;
    end

    function applyIncrement(~, delta, scale)
        x = x + scale * delta;
    end

    function measurement = measure()
        err = abs(1 - x);
        measurement = struct('errorVec', err, 'converged', false, 'scalarError', err);
    end

    function validation = validateCandidate(~, ~, ~)
        validation = struct('isValid', false, 'message', 'always invalid');
    end
end

function testUnknownOptionFieldFailsFast(testCase)
    x = 0;

    problem = struct();
    problem.linearize = @linearize;
    problem.snapshot = @snapshot;
    problem.restore = @restore;
    problem.applyIncrement = @applyIncrement;
    problem.measure = @measure;

    verifyError(testCase, ...
        @() solve_generic_newton(problem, struct('maxIteratoins', 1)), ...
        'solve_generic_newton:UnknownOption');

    function [JCurrent, residual, meta] = linearize()
        JCurrent = 1;
        residual = 1 - x;
        meta = struct('iterateNorm', abs(x));
    end

    function value = snapshot()
        value = x;
    end

    function restore(snapshot)
        x = snapshot;
    end

    function applyIncrement(~, delta, scale)
        x = x + scale * delta;
    end

    function measurement = measure()
        err = abs(1 - x);
        measurement = struct('errorVec', err, 'converged', false, 'scalarError', err);
    end
end
