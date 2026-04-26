function tests = test_solve_regularized_linear_system
%TEST_SOLVE_REGULARIZED_LINEAR_SYSTEM Shared direct/LM step solver coverage.

    tests = functiontests(localfunctions);
end

function testBestEffortUsesDirectOnWellConditionedSystem(testCase)
    A = [2, 0; 0, 1];
    b = [2; 1];

    result = solve_regularized_linear_system(A, b, default_options(), 'best_effort');

    verifyTrue(testCase, result.success);
    verifyEqual(testCase, result.solver, 'direct');
    verifyEqual(testCase, result.lambda, 0, 'AbsTol', 0);
    verifyEqual(testCase, A * result.solution, b, 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, result.conditionEstimate, 1e-10);
    verifyEqual(testCase, result.linearSystemScaling, 'row');
end

function testBestEffortFallsBackToLMOnIllConditionedSystem(testCase)
    A = [1, 0; 1, 1e-12; 0, 1e-12];
    b = [1; 1; 0];
    opts = default_options();
    opts.linearSystemScaling = 'none';

    result = solve_regularized_linear_system(A, b, opts, 'best_effort');

    verifyTrue(testCase, result.success);
    verifyEqual(testCase, result.solver, 'lm');
    verifyGreaterThan(testCase, result.lambda, 0);
    verifyLessThan(testCase, result.directConditionEstimate, opts.directConditionThreshold);
    verifyTrue(testCase, isfinite(result.conditionEstimate));
    verifyEqual(testCase, result.linearSystemScaling, 'none');
end

function testRowScalingCanMakeDirectSolveViable(testCase)
    A = [1e-12, 0; 0, 1];
    b = [1e-12; 2];

    result = solve_regularized_linear_system(A, b, default_options(), 'best_effort');

    verifyTrue(testCase, result.success);
    verifyEqual(testCase, result.solver, 'direct');
    verifyEqual(testCase, result.linearSystemScaling, 'row');
    verifyGreaterThan(testCase, result.directConditionEstimate, ...
        default_options().directConditionThreshold);
    verifyEqual(testCase, A * result.solution, b, 'AbsTol', 1e-12);
end

function testScalingNonePreservesIllConditionedDirectRejection(testCase)
    A = [1e-12, 0; 0, 1];
    b = [1e-12; 2];
    opts = default_options();
    opts.linearSystemScaling = 'none';

    result = solve_regularized_linear_system(A, b, opts, 'direct_only');

    verifyFalse(testCase, result.success);
    verifyEqual(testCase, result.linearSystemScaling, 'none');
    verifySubstring(testCase, result.message, 'ill-conditioned');
    verifyLessThan(testCase, result.directConditionEstimate, opts.directConditionThreshold);
end

function testLmOnlyRejectsNonFiniteInputs(testCase)
    A = [1, NaN; 0, 1];
    b = [1; 0];

    result = solve_regularized_linear_system(A, b, default_options(), 'lm_only');

    verifyFalse(testCase, result.success);
    verifySubstring(testCase, result.message, 'non-finite');
end

function testDirectOnlyReportsIllConditionedRectangularSystem(testCase)
    A = [1, 0; 1, 1e-12; 0, 1e-12];
    b = [1; 1; 0];
    opts = default_options();
    opts.linearSystemScaling = 'none';

    result = solve_regularized_linear_system(A, b, opts, 'direct_only');

    verifyFalse(testCase, result.success);
    verifySubstring(testCase, result.message, 'ill-conditioned');
    verifyLessThan(testCase, result.directConditionEstimate, opts.directConditionThreshold);
end

function testUnknownOptionFieldFailsFast(testCase)
    A = eye(2);
    b = [1; 0];

    verifyError(testCase, ...
        @() solve_regularized_linear_system(A, b, struct('lambdaGorw', 10), 'best_effort'), ...
        'solve_regularized_linear_system:UnknownOption');
end

function testInvalidScalingOptionFailsFast(testCase)
    A = eye(2);
    b = [1; 0];
    opts = default_options();
    opts.linearSystemScaling = 'column';

    verifyError(testCase, ...
        @() solve_regularized_linear_system(A, b, opts, 'best_effort'), ...
        'solve_regularized_linear_system:InvalidLinearSystemScaling');
end

function opts = default_options()
    opts = struct( ...
        'initialLambda', 1e-8, ...
        'lambdaMin', 1e-12, ...
        'lambdaMax', 1e6, ...
        'lambdaGrow', 10, ...
        'directConditionThreshold', 1e-10, ...
        'lmConditionThreshold', 1e-12, ...
        'linearSystemScaling', 'row');
end
