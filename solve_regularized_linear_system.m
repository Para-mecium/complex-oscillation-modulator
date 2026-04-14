function result = solve_regularized_linear_system(A, b, opts, mode)
%SOLVE_REGULARIZED_LINEAR_SYSTEM Shared direct/LM linear solve primitive.

    if nargin < 4 || isempty(mode)
        mode = 'best_effort';
    end

    opts = normalize_options(opts);
    mode = normalize_mode(mode);
    result = empty_result(A);

    [isValid, validationMessage, rhs] = validate_inputs(A, b);
    if ~isValid
        result.message = validationMessage;
        return
    end

    directCondition = estimate_linear_condition(A);
    result.directConditionEstimate = directCondition;

    if ~strcmp(mode, 'lm_only')
        directResult = attempt_direct_solve(A, rhs, directCondition, opts, result);
        if directResult.success || strcmp(mode, 'direct_only')
            result = directResult;
            return
        end
        result.message = directResult.message;
    end

    if strcmp(mode, 'direct_only')
        return
    end

    lmResult = attempt_lm_solve(A, rhs, opts, result);
    result = lmResult;
end

function result = attempt_direct_solve(A, rhs, directCondition, opts, result)
    if ~isfinite(directCondition)
        result.message = 'direct solve condition estimate is not finite';
        return
    end

    if directCondition < opts.directConditionThreshold
        result.message = 'direct solve rejected by ill-conditioned system estimate';
        return
    end

    try
        solution = A \ rhs;
    catch ME
        result.message = format_exception_message('direct solve failed', ME);
        return
    end

    stepNorm = norm(solution, inf);
    if ~all(isfinite(solution)) || ~isfinite(stepNorm)
        result.message = 'direct solve produced a non-finite solution';
        return
    end

    result.success = true;
    result.solution = solution;
    result.solver = 'direct';
    result.conditionEstimate = directCondition;
    result.lambda = 0;
    result.stepNorm = stepNorm;
    result.message = 'direct solve accepted';
end

function result = attempt_lm_solve(A, rhs, opts, result)
    gram = A' * A;
    normalRhs = A' * rhs;
    diagonalScale = normal_equation_scale(gram);
    usedLambda = clamp_lambda(opts.initialLambda, opts);

    while true
        dampingMatrix = gram + usedLambda * diag(diagonalScale);
        conditionEstimate = safe_rcond(dampingMatrix);

        if all(isfinite(dampingMatrix(:))) && all(isfinite(normalRhs)) && ...
                isfinite(conditionEstimate) && conditionEstimate >= opts.lmConditionThreshold
            try
                solution = dampingMatrix \ normalRhs;
            catch ME
                result.message = format_exception_message('LM solve failed', ME);
                return
            end

            stepNorm = norm(solution, inf);
            if all(isfinite(solution)) && isfinite(stepNorm)
                result.success = true;
                result.solution = solution;
                result.solver = 'lm';
                result.conditionEstimate = conditionEstimate;
                result.lambda = usedLambda;
                result.stepNorm = stepNorm;
                result.message = 'LM solve accepted';
                return
            end
        end

        if usedLambda >= opts.lambdaMax
            result.lambda = usedLambda;
            result.message = 'unable to regularize linear system';
            return
        end

        usedLambda = grow_lambda(usedLambda, opts);
    end
end

function value = estimate_linear_condition(A)
    if ~ismatrix(A) || isempty(A) || ~all(isfinite(A(:)))
        value = NaN;
        return
    end

    if size(A, 1) == size(A, 2)
        value = safe_rcond(A);
        return
    end

    gram = A' * A;
    gramCondition = safe_rcond(gram);
    if ~isfinite(gramCondition) || gramCondition < 0
        value = NaN;
    else
        value = sqrt(gramCondition);
    end
end

function scale = normal_equation_scale(gram)
    diagonalScale = diag(gram);
    if isempty(diagonalScale)
        scale = ones(size(gram, 1), 1);
        return
    end

    scale = max(abs(diagonalScale), 1);
end

function [isValid, message, rhs] = validate_inputs(A, b)
    rhs = [];
    isValid = false;
    message = '';

    if ~isnumeric(A) || ~ismatrix(A)
        message = 'coefficient matrix must be a numeric 2-D matrix';
        return
    end

    if ~isnumeric(b)
        message = 'right-hand side must be numeric';
        return
    end

    rhs = reshape(b, [], 1);
    if size(A, 1) ~= numel(rhs)
        message = 'linear system dimensions are inconsistent';
        return
    end

    if isempty(A) || any(~isfinite(A(:))) || any(~isfinite(rhs))
        message = 'linear system contains non-finite entries';
        return
    end

    isValid = true;
end

function mode = normalize_mode(mode)
    mode = char(string(mode));
    validModes = {'best_effort', 'direct_only', 'lm_only'};
    if ~any(strcmp(mode, validModes))
        error('solve_regularized_linear_system:InvalidMode', ...
            'mode must be one of: %s.', strjoin(validModes, ', '));
    end
end

function opts = normalize_options(opts)
    defaults = default_options();
    if nargin == 0 || isempty(opts)
        opts = defaults;
    else
        if ~isstruct(opts)
            error('solve_regularized_linear_system:InvalidOptions', 'opts must be a struct.');
        end
        names = fieldnames(opts);
        for i = 1:numel(names)
            defaults.(names{i}) = opts.(names{i});
        end
        opts = defaults;
    end

    validateattributes(opts.initialLambda, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_regularized_linear_system', 'initialLambda');
    validateattributes(opts.lambdaMin, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_regularized_linear_system', 'lambdaMin');
    validateattributes(opts.lambdaMax, {'numeric'}, {'scalar', 'positive', '>=', opts.lambdaMin}, ...
        'solve_regularized_linear_system', 'lambdaMax');
    validateattributes(opts.lambdaGrow, {'numeric'}, {'scalar', '>', 1}, ...
        'solve_regularized_linear_system', 'lambdaGrow');
    validateattributes(opts.directConditionThreshold, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_regularized_linear_system', 'directConditionThreshold');
    validateattributes(opts.lmConditionThreshold, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_regularized_linear_system', 'lmConditionThreshold');
end

function opts = default_options()
    opts = struct();
    opts.initialLambda = 1e-8;
    opts.lambdaMin = 1e-12;
    opts.lambdaMax = 1e12;
    opts.lambdaGrow = 10;
    opts.directConditionThreshold = 1e-10;
    opts.lmConditionThreshold = 1e-12;
end

function result = empty_result(A)
    result = struct( ...
        'success', false, ...
        'solution', zeros(size(A, 2), 1), ...
        'solver', '', ...
        'conditionEstimate', NaN, ...
        'directConditionEstimate', NaN, ...
        'lambda', NaN, ...
        'stepNorm', NaN, ...
        'message', '');
end

function value = safe_rcond(matrixValue)
    if isempty(matrixValue) || size(matrixValue, 1) ~= size(matrixValue, 2) || ...
            ~all(isfinite(matrixValue(:)))
        value = NaN;
        return
    end

    try
        value = rcond(matrixValue);
    catch
        value = NaN;
    end
end

function value = clamp_lambda(value, opts)
    value = min(max(value, opts.lambdaMin), opts.lambdaMax);
end

function value = grow_lambda(value, opts)
    nextValue = min(opts.lambdaMax, value * opts.lambdaGrow);
    if nextValue <= value
        value = opts.lambdaMax;
    else
        value = nextValue;
    end
end

function message = format_exception_message(prefix, ME)
    if isempty(ME.identifier)
        message = sprintf('%s: %s', prefix, ME.message);
    else
        message = sprintf('%s: %s', prefix, ME.identifier);
    end
end
