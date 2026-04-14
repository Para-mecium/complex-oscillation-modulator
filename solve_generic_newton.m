function result = solve_generic_newton(problem, opts)
%SOLVE_GENERIC_NEWTON Condition-guided direct Newton/LM solver.

    validate_problem(problem);
    opts = normalize_options(opts);

    [J, r, meta] = problem.linearize();
    r = r(:);
    measurement = normalize_measurement(problem.measure());
    [objective, residualNorm] = residual_objective(r);

    history = repmat(empty_history_entry(), 1, opts.maxIterations);
    iterations = 0;
    converged = measurement.converged;
    lambdaSeed = clamp_lambda(opts.initialLambda, opts);

    if converged
        message = 'initial error satisfied tolerance';
    elseif ~all(isfinite(r)) || ~isfinite(objective)
        message = 'initial linearization produced a non-finite residual';
    else
        message = 'maximum iterations reached';
    end

    for iteration = 1:opts.maxIterations
        if converged || ~all(isfinite(r)) || ~isfinite(objective)
            break
        end

        iterationResult = run_iteration(problem, J, r, meta, measurement, lambdaSeed, opts);
        iterations = iteration;
        history(iteration) = make_history_entry(iteration, iterationResult);

        switch iterationResult.status
            case 'accepted'
                J = iterationResult.J;
                r = iterationResult.r;
                meta = iterationResult.meta;
                measurement = iterationResult.measurement;
                objective = iterationResult.objective;
                residualNorm = iterationResult.residualNorm;
                lambdaSeed = iterationResult.nextLambda;

                if measurement.converged
                    converged = true;
                    message = 'error tolerance satisfied';
                    break
                end

            case 'step_too_small'
                message = iterationResult.message;
                lambdaSeed = iterationResult.nextLambda;
                break

            otherwise
                message = iterationResult.message;
                lambdaSeed = iterationResult.nextLambda;
                break
        end
    end

    if iterations == 0
        history = history([]);
    else
        history = history(1:iterations);
    end

    result = struct();
    result.converged = converged;
    result.iterations = iterations;
    result.message = message;
    result.history = history;
    result.errorVec = measurement.errorVec;
    result.scalarError = measurement.scalarError;
    result.linearResidual = r;
    result.linearResidualNorm = residualNorm;
    result.objective = objective;
    result.stepAccepted = ~isempty(history) && any([history.accepted]);
end

function iterationResult = run_iteration(problem, J, r, meta, measurement, lambdaSeed, opts)
    iterationResult = empty_iteration_result(J, r, meta, measurement, lambdaSeed, opts);
    baseSnapshot = problem.snapshot();
    directAttempt = attempt_direct_step(problem, baseSnapshot, J, r, meta, opts, iterationResult);
    if ~strcmp(directAttempt.status, 'failed')
        iterationResult = directAttempt;
        return
    end

    lmAttempt = attempt_lm_step(problem, baseSnapshot, J, r, meta, lambdaSeed, opts, directAttempt);
    if ~strcmp(lmAttempt.status, 'failed')
        iterationResult = lmAttempt;
        return
    end

    iterationResult = lmAttempt;
    iterationResult.message = compose_failure_message(directAttempt.message, lmAttempt.message);
end

function attempt = attempt_direct_step(problem, baseSnapshot, J, r, meta, opts, template)
    attempt = template;
    attempt.solver = 'direct';
    attempt.lambda = 0;
    attempt.nextLambda = clamp_lambda(opts.initialLambda, opts);
    solveResult = solve_regularized_linear_system(J, r, opts, 'direct_only');
    attempt.conditionEstimate = solveResult.conditionEstimate;
    attempt.directConditionEstimate = solveResult.directConditionEstimate;

    if ~solveResult.success
        attempt.message = solveResult.message;
        return
    end

    delta = solveResult.solution;
    attempt.stepNorm = solveResult.stepNorm;

    if attempt.stepNorm <= opts.incrementTolerance
        attempt.status = 'step_too_small';
        attempt.message = 'increment tolerance satisfied';
        return
    end

    [accepted, candidate, candidateMessage, candidateStatus, acceptedScale, backtracks, scaledStepNorm] = ...
        apply_candidate(problem, baseSnapshot, meta, delta, attempt.stepNorm, opts);
    if accepted
        attempt = accept_candidate( ...
            attempt, candidate, opts.initialLambda, opts, acceptedScale, backtracks, scaledStepNorm);
    else
        attempt.status = candidateStatus;
        attempt.message = candidateMessage;
        attempt.acceptedScale = acceptedScale;
        attempt.backtracks = backtracks;
        attempt.stepNorm = scaledStepNorm;
    end
end

function attempt = attempt_lm_step(problem, baseSnapshot, J, r, meta, lambdaSeed, opts, template)
    attempt = template;
    attempt.solver = 'lm';
    lmOpts = opts;
    lmOpts.initialLambda = lambdaSeed;
    solveResult = solve_regularized_linear_system(J, r, lmOpts, 'lm_only');
    attempt.lambda = solveResult.lambda;
    attempt.nextLambda = solveResult.lambda;
    attempt.conditionEstimate = solveResult.conditionEstimate;
    attempt.directConditionEstimate = solveResult.directConditionEstimate;
    attempt.stepNorm = solveResult.stepNorm;

    if ~solveResult.success
        attempt.message = solveResult.message;
        return
    end

    delta = solveResult.solution;

    if attempt.stepNorm <= opts.incrementTolerance
        attempt.status = 'step_too_small';
        attempt.message = 'increment tolerance satisfied';
        return
    end

    [accepted, candidate, candidateMessage, candidateStatus, acceptedScale, backtracks, scaledStepNorm] = ...
        apply_candidate(problem, baseSnapshot, meta, delta, attempt.stepNorm, opts);
    if accepted
        attempt = accept_candidate( ...
            attempt, candidate, max(opts.lambdaMin, solveResult.lambda * opts.lambdaShrink), ...
            opts, acceptedScale, backtracks, scaledStepNorm);
    else
        attempt.status = candidateStatus;
        attempt.message = candidateMessage;
        attempt.acceptedScale = acceptedScale;
        attempt.backtracks = backtracks;
        attempt.stepNorm = scaledStepNorm;
    end
end

function [accepted, candidate, message, status, acceptedScale, backtracks, scaledStepNorm] = ...
        apply_candidate(problem, baseSnapshot, meta, delta, stepNorm, opts)
    accepted = false;
    candidate = struct();
    message = '';
    status = 'failed';
    acceptedScale = 0;
    backtracks = 0;
    scaledStepNorm = stepNorm;

    hasValidator = isfield(problem, 'validateCandidate') && ~isempty(problem.validateCandidate);
    scale = 1;
    while true
        scaledStepNorm = scale * stepNorm;
        if scaledStepNorm <= opts.incrementTolerance
            problem.restore(baseSnapshot);
            message = 'increment tolerance satisfied';
            status = 'step_too_small';
            return
        end

        problem.restore(baseSnapshot);
        try
            problem.applyIncrement(meta, delta, scale);
            if hasValidator
                validation = normalize_candidate_validation(problem.validateCandidate(meta, delta, scale));
                if ~validation.isValid
                    problem.restore(baseSnapshot);
                    message = compose_candidate_validation_message(validation.message);
                    if backtracks >= opts.candidateBacktrackingMaxBacktracks
                        return
                    end
                    backtracks = backtracks + 1;
                    scale = scale * opts.candidateBacktrackingFactor;
                    continue
                end
            end

            [JCandidate, rCandidate, metaCandidate] = problem.linearize();
            rCandidate = rCandidate(:);
            measurementCandidate = normalize_measurement(problem.measure());
            [objectiveCandidate, residualNormCandidate] = residual_objective(rCandidate);
        catch ME
            problem.restore(baseSnapshot);
            message = format_exception_message(ME);
            return
        end

        if ~all(isfinite(rCandidate)) || ~all(isfinite(measurementCandidate.errorVec)) || ...
                ~isfinite(objectiveCandidate) || ~isfinite(measurementCandidate.scalarError)
            problem.restore(baseSnapshot);
            message = 'trial iterate produced non-finite residuals or errors';
            return
        end

        accepted = true;
        status = 'accepted';
        acceptedScale = scale;
        candidate.J = JCandidate;
        candidate.r = rCandidate;
        candidate.meta = metaCandidate;
        candidate.measurement = measurementCandidate;
        candidate.objective = objectiveCandidate;
        candidate.residualNorm = residualNormCandidate;
        return
    end
end

function attempt = accept_candidate(attempt, candidate, nextLambda, opts, acceptedScale, backtracks, scaledStepNorm)
    attempt.status = 'accepted';
    attempt.message = 'accepted';
    attempt.J = candidate.J;
    attempt.r = candidate.r;
    attempt.meta = candidate.meta;
    attempt.measurement = candidate.measurement;
    attempt.objective = candidate.objective;
    attempt.residualNorm = candidate.residualNorm;
    attempt.acceptedScale = acceptedScale;
    attempt.backtracks = backtracks;
    attempt.stepNorm = scaledStepNorm;
    attempt.nextLambda = clamp_lambda(nextLambda, opts);
end


function value = clamp_lambda(value, opts)
    value = min(max(value, opts.lambdaMin), opts.lambdaMax);
end

function [objective, residualNorm] = residual_objective(r)
    residualNorm = norm(r, 2);
    objective = 0.5 * (residualNorm ^ 2);
end

function measurement = normalize_measurement(measurement)
    if ~isstruct(measurement) || ~all(isfield(measurement, {'errorVec', 'converged'}))
        error('solve_generic_newton:InvalidMeasurement', ...
            'problem.measure() must return a struct with fields errorVec and converged.')
    end

    measurement.errorVec = reshape(measurement.errorVec, 1, []);
    if isfield(measurement, 'scalarError') && ~isempty(measurement.scalarError)
        measurement.scalarError = measurement.scalarError;
    elseif isempty(measurement.errorVec)
        measurement.scalarError = 0;
    else
        measurement.scalarError = max(measurement.errorVec);
    end
    measurement.converged = logical(measurement.converged);
end

function entry = make_history_entry(iteration, iterationResult)
    entry = empty_history_entry();
    entry.iteration = iteration;
    entry.objective = iterationResult.objective;
    entry.residualNorm = iterationResult.residualNorm;
    entry.stepNorm = iterationResult.stepNorm;
    entry.acceptedScale = iterationResult.acceptedScale;
    entry.lambda = iterationResult.lambda;
    entry.backtracks = iterationResult.backtracks;
    entry.accepted = strcmp(iterationResult.status, 'accepted');
    entry.maxError = iterationResult.measurement.scalarError;
    entry.solver = iterationResult.solver;
    entry.conditionEstimate = iterationResult.conditionEstimate;
    entry.directConditionEstimate = iterationResult.directConditionEstimate;
end

function entry = empty_history_entry()
    entry = struct( ...
        'iteration', 0, ...
        'objective', NaN, ...
        'residualNorm', NaN, ...
        'stepNorm', NaN, ...
        'acceptedScale', 0, ...
        'lambda', NaN, ...
        'backtracks', 0, ...
        'accepted', false, ...
        'maxError', NaN, ...
        'solver', '', ...
        'conditionEstimate', NaN, ...
        'directConditionEstimate', NaN);
end

function iterationResult = empty_iteration_result(J, r, meta, measurement, lambdaSeed, opts)
    [objective, residualNorm] = residual_objective(r);

    iterationResult = struct();
    iterationResult.status = 'failed';
    iterationResult.message = 'unable to produce a valid Newton update';
    iterationResult.J = J;
    iterationResult.r = r;
    iterationResult.meta = meta;
    iterationResult.measurement = measurement;
    iterationResult.objective = objective;
    iterationResult.residualNorm = residualNorm;
    iterationResult.lambda = clamp_lambda(lambdaSeed, opts);
    iterationResult.nextLambda = clamp_lambda(lambdaSeed, opts);
    iterationResult.acceptedScale = 0;
    iterationResult.backtracks = 0;
    iterationResult.stepNorm = 0;
    iterationResult.solver = '';
    iterationResult.conditionEstimate = NaN;
    iterationResult.directConditionEstimate = NaN;
end

function message = compose_failure_message(primaryMessage, secondaryMessage)
    if isempty(primaryMessage)
        message = secondaryMessage;
        return
    end
    if isempty(secondaryMessage)
        message = primaryMessage;
        return
    end

    message = sprintf('%s; fallback LM failed: %s', primaryMessage, secondaryMessage);
end

function message = format_exception_message(ME)
    if ~isempty(ME.identifier)
        message = sprintf('trial iterate rejected: %s', ME.identifier);
    else
        message = sprintf('trial iterate rejected: %s', ME.message);
    end
end

function validation = normalize_candidate_validation(validation)
    if ~isstruct(validation) || ~isfield(validation, 'isValid')
        error('solve_generic_newton:InvalidCandidateValidation', ...
            ['problem.validateCandidate() must return a struct with field ' ...
            '''isValid''.'])
    end

    validation.isValid = logical(validation.isValid);
    if ~isfield(validation, 'message') || isempty(validation.message)
        validation.message = '';
    elseif isstring(validation.message) && isscalar(validation.message)
        validation.message = char(validation.message);
    end
end

function message = compose_candidate_validation_message(message)
    if isempty(message)
        message = 'candidate validation failed';
        return
    end

    message = sprintf('candidate validation failed: %s', message);
end

function validate_problem(problem)
    requiredFields = {'linearize', 'snapshot', 'restore', 'applyIncrement', 'measure'};
    if ~isstruct(problem) || ~all(isfield(problem, requiredFields))
        error('solve_generic_newton:InvalidProblem', ...
            'problem must be a struct with fields linearize, snapshot, restore, applyIncrement, and measure.')
    end

    if isfield(problem, 'validateCandidate') && ~isempty(problem.validateCandidate) && ...
            ~isa(problem.validateCandidate, 'function_handle')
        error('solve_generic_newton:InvalidProblem', ...
            'problem.validateCandidate must be a function handle when provided.')
    end
end

function opts = normalize_options(opts)
    defaults = default_options();
    if nargin == 0 || isempty(opts)
        opts = defaults;
    else
        if ~isstruct(opts)
            error('solve_generic_newton:InvalidOptions', 'opts must be a struct.')
        end
        names = fieldnames(opts);
        for i = 1:numel(names)
            defaults.(names{i}) = opts.(names{i});
        end
        opts = defaults;
    end

    validateattributes(opts.maxIterations, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
        'solve_generic_newton', 'maxIterations');
    validateattributes(opts.incrementTolerance, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_generic_newton', 'incrementTolerance');
    validateattributes(opts.initialLambda, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_generic_newton', 'initialLambda');
    validateattributes(opts.lambdaMin, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_generic_newton', 'lambdaMin');
    validateattributes(opts.lambdaMax, {'numeric'}, {'scalar', 'positive', '>=', opts.lambdaMin}, ...
        'solve_generic_newton', 'lambdaMax');
    validateattributes(opts.lambdaGrow, {'numeric'}, {'scalar', '>', 1}, ...
        'solve_generic_newton', 'lambdaGrow');
    validateattributes(opts.lambdaShrink, {'numeric'}, {'scalar', 'positive', '<=', 1}, ...
        'solve_generic_newton', 'lambdaShrink');
    validateattributes(opts.directConditionThreshold, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_generic_newton', 'directConditionThreshold');
    validateattributes(opts.lmConditionThreshold, {'numeric'}, {'scalar', 'positive'}, ...
        'solve_generic_newton', 'lmConditionThreshold');
    validateattributes(opts.candidateBacktrackingFactor, {'numeric'}, ...
        {'scalar', 'positive', '<', 1}, 'solve_generic_newton', 'candidateBacktrackingFactor');
    validateattributes(opts.candidateBacktrackingMaxBacktracks, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'}, ...
        'solve_generic_newton', 'candidateBacktrackingMaxBacktracks');
end

function opts = default_options()
    opts = struct();
    opts.maxIterations = 18;
    opts.incrementTolerance = 1e-8;
    opts.initialLambda = 1e-8;
    opts.lambdaMin = 1e-12;
    opts.lambdaMax = 1e12;
    opts.lambdaGrow = 10;
    opts.lambdaShrink = 0.3;
    opts.directConditionThreshold = 1e-10;
    opts.lmConditionThreshold = 1e-12;
    opts.candidateBacktrackingFactor = 0.5;
    opts.candidateBacktrackingMaxBacktracks = 6;
end
