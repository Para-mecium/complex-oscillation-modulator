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
end

function testBestEffortFallsBackToLMOnIllConditionedSystem(testCase)
    A = [1, 0; 1, 1e-12; 0, 1e-12];
    b = [1; 1; 0];

    result = solve_regularized_linear_system(A, b, default_options(), 'best_effort');

    verifyTrue(testCase, result.success);
    verifyEqual(testCase, result.solver, 'lm');
    verifyGreaterThan(testCase, result.lambda, 0);
    verifyLessThan(testCase, result.directConditionEstimate, default_options().directConditionThreshold);
    verifyTrue(testCase, isfinite(result.conditionEstimate));
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

    result = solve_regularized_linear_system(A, b, default_options(), 'direct_only');

    verifyFalse(testCase, result.success);
    verifySubstring(testCase, result.message, 'ill-conditioned');
    verifyLessThan(testCase, result.directConditionEstimate, default_options().directConditionThreshold);
end

function opts = default_options()
    opts = struct( ...
        'initialLambda', 1e-8, ...
        'lambdaMin', 1e-12, ...
        'lambdaMax', 1e6, ...
        'lambdaGrow', 10, ...
        'directConditionThreshold', 1e-10, ...
        'lmConditionThreshold', 1e-12);
end
