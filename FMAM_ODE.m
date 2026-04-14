classdef FMAM_ODE < handle
    properties (SetAccess = private)
        sys % the functions of the dynamical system, expected to be a 1*N cell
        obs % the functions of observables, expected to be a 1*n cell
        derivatives % the derivatives of system function and observables
        items_perturb % quantities which need to be modulated, expected to be a 1*n_per struct 
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %                 Attributes                   %
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %   -- prop:      propname(str)    
                      %
                      %      supported proplist = {params, p_Psi, q_Psi,
                      %                  p_var, q_var, varPhiMax,
                      %                  obsPhiMax, varPhiMin, obsPhiMin,
                      %                  varAmp, obsAmp, varMax, obsMax,
                      %                  varMin, obsMin, varPhase}
                      %
                      %   -- idx:       index(int)
                      %
                      %      format:    p_var(q_var)
                      %                      -- [i,j] for jth the Fourier
                      %                      coefficient of ith term
                      %                 varPhase
                      %                      -- [i,j] for the phase difference 
                      %                         between var_i(obs_i) and 
                      %                         var_j(obs_j)
                      %                 others -- i 
                      %
                      %   -- target:    modulation target(double)
                      %   For props that have only one single component,
                      %   the 'idx' attr should be set as 1.
                      %
                      %   E.g. iterms_perturb(1,i) = struct('prop',
                      %   'varAmp', 'idx', 4, 'target', 10) means that the
                      %   ith modulation target is to vary the amplitude of
                      %   the 4th system variable to 10.

        items_controlled % quantities which is varied to achieve the modulation target, 
                         % expected to be a vector
        maxstepsize % ith component is the maxstepsize of ith quantity of items_perturb
        truncationOrder % truncation order of trigonometric polynomials
        accuracy % 
        newtonOptions = struct()
        continuationOptions = struct()
        logs % used to save the info of each steps
        continuationStatus = struct('completed', false, 'lambda', 0, ...
            'reason', '', 'triggerValue', NaN)
    end

    properties
        errBound = 1e-8;
        Lconst = 500;
        needLog = false
        verbose = true
        checkPsiNonnegative = true
    end

    properties (Access = private)
        dimVar
        dimObs
        dimParams 
        n_per
        psiUpdateMode = false
        targetCurr = []
        stepsize = 0
        targetPath = struct([])
        
        p_Psi_init = []
        q_Psi_init = []
        solverView = struct()
    end

    properties (Dependent)
        isPsiUpdated
        items_per_curr
    end

    methods
        function obj = FMAM_ODE(system,observables,initialSolverView,items_per,Coe_Controlled,maxstepsize,err,varargin)
            ctorOptions = FMAM_ODE.parseConstructorOptions(varargin{:});
            if iscell(system)
                system = reshape(system,1,[]);
            end
            if isempty(observables)
                observables = {};
            elseif iscell(observables)
                observables = reshape(observables,1,[]);
            end
            obj.sys = system;
            obj.obs = observables;
            obj.dimVar = numel(system);
            obj.dimObs = numel(observables);

            candidateView = obj.coerceSolverViewInput(initialSolverView);
            obj.truncationOrder = candidateView.truncationOrder;
            obj.dimParams = numel(candidateView.params);
            obj.validateStateDimensions(candidateView);
            obj.solverView = obj.normalizeSolverView(candidateView);

            obj.n_per = numel(items_per);
            [items_per,Coe_Controlled,maxstepsize,err] = obj.validateModulationSetup(items_per,Coe_Controlled,maxstepsize,err);
            obj.items_perturb = items_per;
            obj.items_controlled = Coe_Controlled;
            obj.maxstepsize = maxstepsize;
            obj.accuracy = err;

            obj.derivatives = ctorOptions.derivatives;
            obj.newtonOptions = obj.normalizeNewtonOptions(ctorOptions.newtonOptions, 'newtonOptions');
            obj.continuationOptions = ctorOptions.continuationOptions;
            obj.setPsiUpdateMode(ctorOptions.isPsiUpdated);
            obj.needLog = ctorOptions.needLog;
            obj.verbose = ctorOptions.verbose;
            obj.errBound = ctorOptions.errBound;
            obj.Lconst = ctorOptions.Lconst;
            obj.checkPsiNonnegative = ctorOptions.checkPsiNonnegative;

            obj.continuationOptions = obj.normalizeContinuationOptions(obj.continuationOptions);
            obj.validateDerivativeCache();

            obj.targetCurr = obj.initialTargetValues();
        end

        function val = get.isPsiUpdated(obj)
            val = obj.psiUpdateMode;
        end

        function set.isPsiUpdated(obj,value)
            obj.setPsiUpdateMode(value);
        end

        function val = get.items_per_curr(obj)
            val = obj.targetCurr;
        end

        function setPsiUpdateMode(obj,value)
            obj.applyPsiUpdateMode(value);
        end

        function perturb(obj)
            obj.targetCurr = obj.targetCurr + obj.stepsize;
        end

        function rebuilt = rebuildState(obj)
            rebuilt = state.fromViews(obj.obs,obj.solverView,obj.exportDerivedView());
            rebuilt.checkPsiNonnegative = obj.checkPsiNonnegative;
        end

        function step(obj)
            obj.targetPath = obj.buildContinuationPath();
            obj.continuationStatus = struct( ...
                'completed', false, ...
                'lambda', 0, ...
                'reason', '', ...
                'triggerValue', NaN);
            obj.logMessage('STEP', ...
                'start continuation: nTargets=%d, predictorMode=%s', ...
                obj.n_per, obj.continuationOptions.predictorMode);

            if obj.n_per == 0 || all(abs([obj.targetPath.deltaValue]) <= obj.accuracy)
                obj.stepsize = zeros(1,obj.n_per);
                obj.targetCurr = obj.initialTargetValues();
                obj.continuationStatus = struct( ...
                    'completed', true, ...
                    'lambda', 1, ...
                    'reason', 'no_active_target', ...
                    'triggerValue', NaN);
                obj.logMessage('STEP', ...
                    'no active target (all delta <= accuracy); skip continuation.');
                return
            end

            lambda = 0;
            acceptedState = obj.captureAcceptedContinuationState(lambda,0,struct(),struct());
            acceptedHistory = acceptedState.point;

            dlambdaCap = obj.computeLambdaStepCap(obj.targetPath);
            dlambda = obj.chooseInitialLambdaStep(dlambdaCap);
            obj.logMessage('STEP', ...
                'initial dlambda=%.6e (cap=%.6e, min=%.6e, max=%.6e)', ...
                dlambda, dlambdaCap, obj.continuationOptions.minLambdaStep, ...
                obj.continuationOptions.maxLambdaStep);
            count = 0;
            failCount = 0;
            baseFailureCount = 0;
            failureDiagnostics = struct( ...
                'trialIteration', {}, ...
                'lambdaTrial', {}, ...
                'dlambdaTrial', {}, ...
                'fitResult', {});

            while lambda < 1
                if dlambda < obj.continuationOptions.minLambdaStep
                    obj.logMessage('STOP', ...
                        ['dlambda=%.6e < minLambdaStep=%.6e, terminate continuation at ' ...
                        'lambda=%.6f'], ...
                        dlambda, obj.continuationOptions.minLambdaStep, lambda);
                    break
                end

                dlambdaTrial = obj.clampLambdaStep(dlambda,dlambdaCap,1 - lambda);
                if dlambdaTrial < obj.continuationOptions.minLambdaStep
                    obj.logMessage('STOP', ...
                        ['trial dlambda=%.6e < minLambdaStep=%.6e after clamp, terminate ' ...
                        'continuation at lambda=%.6f'], ...
                        dlambdaTrial, obj.continuationOptions.minLambdaStep, lambda);
                    break
                end

                lambdaTrial = min(1,lambda + dlambdaTrial);
                trialTarget = obj.evaluateTargetPath(obj.targetPath,lambdaTrial);

                predictorInfo = obj.preparePredictorTrial( ...
                    acceptedState,acceptedHistory,lambdaTrial,dlambdaTrial,baseFailureCount);
                obj.activateTrialContinuationState( ...
                    acceptedState,trialTarget,predictorInfo.snapshot);

                count = count + 1;
                obj.logMessage('TRIAL', ...
                    ['iter=%d, lambda %.6f -> %.6f (dlambda=%.6e), predictor=%s, ' ...
                    'fallback=%d, info=%s'], ...
                    count, lambda, lambdaTrial, dlambdaTrial, predictorInfo.used, ...
                    predictorInfo.fallbackCount, predictorInfo.message);

                try
                    fitResult = obj.fit();
                catch ME
                    obj.restoreAcceptedContinuationState(acceptedState);
                    obj.logMessage('ERROR', ...
                        'iter=%d raised %s, rollback to lambda=%.6f', ...
                        count, FMAM_ODE.exceptionSummary(ME), lambda);
                    rethrow(ME)
                end

                if fitResult.converged
                    directCondition = obj.extractDirectConditionEstimate(fitResult);
                    psiIssue = obj.convergedStatePsiIssue();
                    if ~psiIssue.isValid
                        [failureDiagnostics, failCount, baseFailureCount, dlambda, stopNow] = ...
                            obj.rejectContinuationTrial( ...
                                failureDiagnostics, failCount, baseFailureCount, ...
                                acceptedState, count, lambda, lambdaTrial, ...
                                dlambdaTrial, dlambdaCap, fitResult, ...
                                sprintf('rejected: converged state has invalid Psi (%s)', ...
                                psiIssue.message));
                        if stopNow
                            break
                        end
                        continue
                    end

                    if obj.shouldStopForCondition(directCondition)
                        obj.restoreAcceptedContinuationState(acceptedState);
                        obj.stepsize = zeros(1,obj.n_per);
                        obj.continuationStatus = struct( ...
                            'completed', false, ...
                            'lambda', lambda, ...
                            'reason', 'condition_stop', ...
                            'triggerValue', directCondition);
                        obj.logMessage('STOP', ...
                            ['iter=%d accepted predictor-corrector point was discarded ' ...
                            'because directCondition=%s < conditionStopRcond=%s ' ...
                            'at lambdaTrial=%.6f'], ...
                            count, obj.formatScalar(directCondition), ...
                            obj.formatScalar(obj.continuationOptions.conditionStopRcond), ...
                            lambdaTrial);
                        break
                    end

                    lambda = lambdaTrial;
                    acceptedState = obj.captureAcceptedContinuationState( ...
                        lambda,dlambdaTrial,fitResult,predictorInfo);
                    acceptedHistory(end + 1) = acceptedState.point; %#ok<AGROW>

                    failCount = 0;
                    baseFailureCount = 0;

                    if obj.needLog
                        obj.appendAcceptedLog(acceptedState.point,fitResult);
                    end

                    if lambda >= 1
                        acceptedState.stepsize = zeros(1,obj.n_per);
                        obj.stepsize = acceptedState.stepsize;
                        obj.continuationStatus = struct( ...
                            'completed', true, ...
                            'lambda', 1, ...
                            'reason', 'target_reached', ...
                            'triggerValue', directCondition);
                        obj.logMessage('DONE', ...
                            ['iter=%d accepted; reached lambda=1.000000. final target ' ...
                            'state committed.'], count);
                        break
                    end

                    nextDlambda = obj.chooseAcceptedLambdaStep(dlambdaTrial,fitResult,dlambdaCap);
                    obj.logMessage('ACCEPT', ...
                        'iter=%d accepted: lambda=%.6f, next dlambda=%.6e', ...
                        count, lambda, nextDlambda);
                    dlambda = nextDlambda;
                    continue
                end

                [failureDiagnostics, failCount, baseFailureCount, dlambda, stopNow] = ...
                    obj.rejectContinuationTrial( ...
                        failureDiagnostics, failCount, baseFailureCount, ...
                        acceptedState, count, lambda, lambdaTrial, ...
                        dlambdaTrial, dlambdaCap, fitResult, fitResult.message);
                if stopNow
                    break
                end
            end

            obj.restoreAcceptedContinuationState(acceptedState);

            continuationConverged = lambda >= 1 - 10 * eps;
            if continuationConverged
                if ~obj.continuationStatus.completed
                    obj.continuationStatus = struct( ...
                        'completed', true, ...
                        'lambda', lambda, ...
                        'reason', 'target_reached', ...
                        'triggerValue', NaN);
                end
                obj.logMessage('STEP', ...
                    'continuation finished at lambda=%.6f with %d accepted points.', ...
                    lambda, max(0,numel(acceptedHistory)-1));
                return
            end

            if isempty(obj.continuationStatus.reason)
                obj.continuationStatus = struct( ...
                    'completed', false, ...
                    'lambda', lambda, ...
                    'reason', 'stopped_early', ...
                    'triggerValue', NaN);
            end
            obj.logMessage('STEP', ...
                ['continuation stopped early at lambda=%.6f with %d accepted points; ' ...
                'emit Newton diagnostics for rejected trials.'], ...
                lambda, max(0,numel(acceptedHistory)-1));
            for i = 1:numel(failureDiagnostics)
                diagnostic = failureDiagnostics(i);
                obj.logMessage('NEWTON', ...
                    'rejected trial iter=%d, lambdaTrial=%.6f, dlambda=%.6e', ...
                    diagnostic.trialIteration, diagnostic.lambdaTrial, diagnostic.dlambdaTrial);
                obj.logNewtonSummary(diagnostic.fitResult);
            end
        end

        function autostepsize(obj)
            obj.targetPath = obj.buildContinuationPath();
            if obj.n_per == 0
                obj.stepsize = zeros(1,0);
                return
            end

            dlambdaCap = obj.computeLambdaStepCap(obj.targetPath);
            dlambda = obj.chooseInitialLambdaStep(dlambdaCap);
            if dlambda <= 0
                obj.stepsize = zeros(1,obj.n_per);
                return
            end

            nextTarget = obj.evaluateTargetPath(obj.targetPath,min(1,dlambda));
            obj.stepsize = nextTarget - obj.targetCurr;
        end

        function result = fit(obj)
            solveResult = obj.runNewtonSolve(struct());
            result = obj.translateNewtonResult(solveResult,'finalError');
        end

        function view = exportSolverView(obj)
            view = obj.solverView;
        end

        function derived = exportDerivedView(obj)
            derived = obj.buildDerivedView(obj.solverView);
        end

        function loadSolverView(obj,view)
            obj.solverView = obj.normalizeSolverView(view);
        end

        function val = res(obj)
            runtimeContext = obj.currentRuntimeContext();
            system = obj.sys;
            solverView = runtimeContext.solverView;
            derived = runtimeContext.derived;
            targetCtx = runtimeContext.targetCtx;

            p_var = solverView.p_var;
            q_var = solverView.q_var;
            p_Psi = solverView.p_Psi;
            q_Psi = solverView.q_Psi;
            parameters = solverView.params;

            L = obj.Lconst;
            N = obj.dimVar;
            % n = obj.dimObs;
            MVar = obj.truncationOrder;
            MPsi = size(q_Psi,1);

            phi = (0:L-1)'*2*pi/L;
            [vcVar,vsVar] = obj.Vec_CS(phi,MVar,L);
            [vcPsi,vsPsi] = obj.Vec_CS(phi,MPsi,L);
            TS_var = vcVar*p_var + vsVar*q_var;
            Psi = vcPsi*p_Psi + vsPsi*q_Psi;
            TS_Dvar = -(vsVar.*(1:MVar))*p_var(2:end,:) + (vcVar(:,2:end).*(1:MVar))*q_var;

            res_sys = 0;
            res_var_phi = 0;
            res_obs_phi = 0;
            res_target = 0;
            for i = 1:N
                residue = FMAM_ODE.residue_system(system,Psi,TS_var,TS_Dvar,parameters,i);
                res_temp = fft(residue)/L;
                res_sys = max(res_sys,max(abs(res_temp(1:MVar+1))));

                [needMax, needMin] = obj.targetNeedsVariableExtrema(i,targetCtx);

                if needMax
                    p = p_var(:,i);
                    q = q_var(:,i);
                    phiExt = solverView.varPhiMax(i);

                    res_var_phi = max(res_var_phi,abs(FMAM_ODE.residue_phi_var(p,q,phiExt)));
                end

                if needMin
                    p = p_var(:,i);
                    q = q_var(:,i);
                    phiExt = solverView.varPhiMin(i);

                    res_var_phi = max(res_var_phi,abs(FMAM_ODE.residue_phi_var(p,q,phiExt)));
                end
            end

            for k = 1:obj.dimObs
                [needMax, needMin] = obj.targetNeedsObservableExtrema(k,targetCtx);

                if needMax
                    phiExt = solverView.obsPhiMax(k);
                    res_obs_phi = max(res_obs_phi,abs(FMAM_ODE.residue_phi_obs(obj.derivatives.obs,p_var,q_var,phiExt,k)));
                end

                if needMin
                    phiExt = solverView.obsPhiMin(k);
                    res_obs_phi = max(res_obs_phi,abs(FMAM_ODE.residue_phi_obs(obj.derivatives.obs,p_var,q_var,phiExt,k)));
                end
            end

            for j = 1:obj.n_per
                item_per = obj.items_perturb(j);
                val_target = obj.targetCurr(j);
                res_temp = val_target - obj.currentTargetValue(item_per,solverView,derived,targetCtx);
                res_target = max(res_target,abs(res_temp));
            end

            val = [res_sys,res_var_phi,res_obs_phi,res_target];
        end

        function result = oneIter(obj)
            result = obj.translateNewtonResult( ...
                obj.runNewtonSolve(struct('maxIterations',1)),'errorVec');
        end    

    end

    methods (Access = private)
        function result = runNewtonSolve(obj,overrideOptions)
            [problem,options] = obj.buildNewtonProblem(overrideOptions);
            result = solve_generic_newton(problem,options);
        end

        function [problem,options] = buildNewtonProblem(obj,overrideOptions)
            if nargin < 2 || isempty(overrideOptions)
                overrideOptions = struct();
            end

            obj.ensurePsiReferenceInitialized();

            problem = struct();
            problem.linearize = @() obj.linearizeNewtonProblem();
            problem.snapshot = @() obj.snapshotSolverView();
            problem.restore = @(snapshot) obj.restoreSolverView(snapshot);
            problem.applyIncrement = @(meta,increments,scale) obj.applyUnknownIncrement(meta,increments,scale);
            problem.measure = @() obj.measureNewtonState();
            problem.validateCandidate = @(meta,delta,scale) obj.validateNewtonCandidate(meta,delta,scale);
            options = obj.effectiveNewtonOptions(overrideOptions);
        end

        function ensurePsiReferenceInitialized(obj)
            if ~obj.isPsiUpdated && (isempty(obj.p_Psi_init) || isempty(obj.q_Psi_init))
                obj.captureFrozenPsiReference();
            end
        end

        function [A,residual,meta] = linearizeNewtonProblem(obj)
            ctx = obj.buildNewtonAssemblyContext();
            [A,residual,indexMap,unknowns] = assemble_newton_linear_system(ctx);

            meta = struct();
            meta.indexMap = indexMap;
            meta.unknowns = unknowns;
            meta.iterateNorm = obj.currentUnknownInfNorm(unknowns);
        end

        function measurement = measureNewtonState(obj)
            err = reshape(obj.res(),1,[]);
            measurement = struct();
            measurement.errorVec = err;
            measurement.converged = all(err <= obj.errBound);
            measurement.scalarError = max(err);
        end

        function validation = validateNewtonCandidate(obj,~,~,~)
            validation = struct('isValid', true, 'message', '');
            if ~obj.checkPsiNonnegative
                return
            end

            issue = fmam_state_ops.positivePsiIssue(obj.currentPsiValues());
            validation.isValid = issue.isValid;
            validation.message = issue.message;
        end

        function opts = defaultNewtonOptions(~)
            opts = struct();
            opts.maxIterations = 50;
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

        function opts = effectiveNewtonOptions(obj,overrideOptions)
            opts = obj.normalizeNewtonOptions(obj.newtonOptions, 'newtonOptions');
            opts = FMAM_ODE.mergeOptionStruct( ...
                opts,overrideOptions,'overrideOptions','FMAM_ODE:UnknownNewtonOption');
            opts = obj.normalizeNewtonOptions(opts, 'overrideOptions');
        end

        function opts = linearSolverOptions(obj,overrideOptions)
            newtonOpts = obj.effectiveNewtonOptions(overrideOptions);
            opts = struct( ...
                'initialLambda', newtonOpts.initialLambda, ...
                'lambdaMin', newtonOpts.lambdaMin, ...
                'lambdaMax', newtonOpts.lambdaMax, ...
                'lambdaGrow', newtonOpts.lambdaGrow, ...
                'directConditionThreshold', newtonOpts.directConditionThreshold, ...
                'lmConditionThreshold', newtonOpts.lmConditionThreshold);
        end

        function opts = normalizeNewtonOptions(obj,options,argName)
            if nargin < 3 || isempty(argName)
                argName = 'newtonOptions';
            end

            opts = obj.defaultNewtonOptions();
            if nargin < 2 || isempty(options)
                return
            end

            opts = FMAM_ODE.mergeOptionStruct( ...
                opts,options,argName,'FMAM_ODE:UnknownNewtonOption');

            validateattributes(opts.maxIterations, {'numeric'}, ...
                {'scalar', 'integer', 'nonnegative'}, 'FMAM_ODE', [argName '.maxIterations']);
            validateattributes(opts.incrementTolerance, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', [argName '.incrementTolerance']);
            validateattributes(opts.initialLambda, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', [argName '.initialLambda']);
            validateattributes(opts.lambdaMin, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', [argName '.lambdaMin']);
            validateattributes(opts.lambdaMax, {'numeric'}, ...
                {'scalar', 'positive', '>=', opts.lambdaMin}, 'FMAM_ODE', [argName '.lambdaMax']);
            validateattributes(opts.lambdaGrow, {'numeric'}, ...
                {'scalar', '>', 1}, 'FMAM_ODE', [argName '.lambdaGrow']);
            validateattributes(opts.lambdaShrink, {'numeric'}, ...
                {'scalar', 'positive', '<=', 1}, 'FMAM_ODE', [argName '.lambdaShrink']);
            validateattributes(opts.directConditionThreshold, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', [argName '.directConditionThreshold']);
            validateattributes(opts.lmConditionThreshold, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', [argName '.lmConditionThreshold']);
            validateattributes(opts.candidateBacktrackingFactor, {'numeric'}, ...
                {'scalar', 'positive', '<', 1}, 'FMAM_ODE', [argName '.candidateBacktrackingFactor']);
            validateattributes(opts.candidateBacktrackingMaxBacktracks, {'numeric'}, ...
                {'scalar', 'integer', 'nonnegative'}, ...
                'FMAM_ODE', [argName '.candidateBacktrackingMaxBacktracks']);
        end

        function ctx = buildNewtonAssemblyContext(obj)
            runtimeContext = obj.currentRuntimeContext();
            solverView = runtimeContext.solverView;
            ctx = struct();
            ctx.sys = obj.sys;
            ctx.obs = obj.obs;
            ctx.derivatives = obj.derivatives;

            ctx.L = obj.Lconst;
            ctx.dimVar = obj.dimVar;
            ctx.dimObs = obj.dimObs;
            ctx.truncationOrder = obj.truncationOrder;
            ctx.dimParams = obj.dimParams;

            ctx.params = solverView.params;
            ctx.p_Psi = solverView.p_Psi;
            ctx.q_Psi = solverView.q_Psi;
            ctx.p_var = solverView.p_var;
            ctx.q_var = solverView.q_var;
            ctx.varPhiMax = solverView.varPhiMax;
            ctx.varPhiMin = solverView.varPhiMin;
            ctx.obsPhiMax = solverView.obsPhiMax;
            ctx.obsPhiMin = solverView.obsPhiMin;
            ctx.PV = solverView.PV;

            ctx.items_perturb = obj.items_perturb;
            ctx.items_controlled = obj.items_controlled;
            ctx.targetCurr = obj.targetCurr;
            ctx.isPsiUpdated = obj.isPsiUpdated;
            ctx.p_Psi_init = obj.p_Psi_init;
            ctx.q_Psi_init = obj.q_Psi_init;
            ctx.targetRuleCtx = runtimeContext.targetCtx;
        end

        function snapshot = snapshotSolverView(obj)
            snapshot = obj.solverView;
        end

        function restoreSolverView(obj,snapshot)
            obj.solverView = obj.normalizeSolverView(snapshot);
        end

        function applyUnknownIncrement(obj,meta,increments,scale)
            if nargin < 4 || isempty(scale)
                scale = 1;
            end
            view = obj.solverView;
            obj.solverView = FMAM_ODE.applyIncrementToSolverView(view,meta,increments,scale);
        end

        function value = currentUnknownInfNorm(obj,unknowns)
            value = 0;
            for i = 1:numel(unknowns)
                currentValue = obj.solverView.(unknowns{i});
                if ~isempty(currentValue)
                    value = max(value,norm(currentValue(:),inf));
                end
            end
        end

        function runtimeContext = currentRuntimeContext(obj,solverView,derived)
            if nargin < 2 || isempty(solverView)
                solverView = obj.solverView;
            end
            if nargin < 3 || isempty(derived)
                derived = obj.buildDerivedView(solverView);
            end
            runtimeContext = struct( ...
                'solverView', solverView, ...
                'derived', derived, ...
                'targetCtx', obj.targetRuleContext(solverView,derived));
        end

        function psiValues = currentPsiValues(obj)
            psiValues = fmam_state_ops.evaluateTrigSeries( ...
                obj.phaseGrid(),obj.solverView.p_Psi,obj.solverView.q_Psi);
        end

        function phi = phaseGrid(~)
            phi = (0:fmam_state_defaults.LphiConst-1)'*2*pi/fmam_state_defaults.LphiConst;
        end

        function solverView = coerceSolverViewInput(~,solverViewInput)
            if ~isstruct(solverViewInput)
                error('FMAM_ODE:InvalidInitialSolverView', ...
                    ['initialSolverView must be a canonical solverView struct. ' ...
                    'Convert state objects with fmam_state_ops.solverViewFromState before constructing FMAM_ODE.'])
            end
            solverView = fmam_state_ops.normalizeSolverView(solverViewInput);
        end

        function normalized = normalizeSolverView(obj,view)
            normalized = fmam_state_ops.normalizeSolverView( ...
                view,obj.dimVar,obj.dimObs,obj.dimParams,obj.truncationOrder);
        end

        function derived = buildDerivedView(obj,solverView)
            derived = fmam_state_ops.buildDerivedView( ...
                obj.obs,solverView, ...
                fmam_state_defaults.LphiConst,fmam_state_defaults.Lconst, ...
                fmam_state_defaults.countMax,fmam_state_defaults.errMax);
        end

        function opts = normalizeContinuationOptions(~,options)
            opts = FMAM_ODE.defaultContinuationOptions();
            if nargin < 2 || isempty(options)
                return
            end
            if ~isstruct(options)
                error('FMAM_ODE:InvalidContinuationOptions', ...
                    'continuationOptions must be a struct.')
            end

            opts = FMAM_ODE.mergeOptionStruct( ...
                opts,options,'continuationOptions','FMAM_ODE:UnknownContinuationOption');

            validateattributes(opts.initialSteps, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}, 'FMAM_ODE', 'continuationOptions.initialSteps');
            if ~(isscalar(opts.initialLambdaStep) && (isnan(opts.initialLambdaStep) || ...
                    (isfinite(opts.initialLambdaStep) && opts.initialLambdaStep > 0)))
                error('FMAM_ODE:InvalidContinuationOptions', ...
                    'continuationOptions.initialLambdaStep must be NaN or a positive finite scalar.')
            end
            validateattributes(opts.minLambdaStep, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', 'continuationOptions.minLambdaStep');
            validateattributes(opts.maxLambdaStep, {'numeric'}, ...
                {'scalar', 'positive', '>=', opts.minLambdaStep}, 'FMAM_ODE', 'continuationOptions.maxLambdaStep');
            validateattributes(opts.growthFactor, {'numeric'}, ...
                {'scalar', '>=', 1}, 'FMAM_ODE', 'continuationOptions.growthFactor');
            validateattributes(opts.shrinkFactor, {'numeric'}, ...
                {'scalar', 'positive', '<=', 1}, 'FMAM_ODE', 'continuationOptions.shrinkFactor');
            validateattributes(opts.backtrackingFactor, {'numeric'}, ...
                {'scalar', 'positive', '<', 1}, 'FMAM_ODE', 'continuationOptions.backtrackingFactor');
            validateattributes(opts.easyIterations, {'numeric'}, ...
                {'scalar', 'integer', 'nonnegative'}, 'FMAM_ODE', 'continuationOptions.easyIterations');
            validateattributes(opts.slowIterations, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}, 'FMAM_ODE', 'continuationOptions.slowIterations');
            validateattributes(opts.maxFailures, {'numeric'}, ...
                {'scalar', 'integer', 'nonnegative'}, 'FMAM_ODE', 'continuationOptions.maxFailures');
            validateattributes(opts.conditionStopEnabled, {'numeric', 'logical'}, ...
                {'scalar'}, 'FMAM_ODE', 'continuationOptions.conditionStopEnabled');
            validateattributes(opts.conditionStopRcond, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', 'continuationOptions.conditionStopRcond');
            validateattributes(opts.quadraticCurvatureThreshold, {'numeric'}, ...
                {'scalar', 'nonnegative'}, 'FMAM_ODE', 'continuationOptions.quadraticCurvatureThreshold');
            validateattributes(opts.hermiteMaxExtrapolationRatio, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', 'continuationOptions.hermiteMaxExtrapolationRatio');
            validateattributes(opts.predictorStepGrowthLimit, {'numeric'}, ...
                {'scalar', 'positive'}, 'FMAM_ODE', 'continuationOptions.predictorStepGrowthLimit');
            validateattributes(opts.quadraticStepRatioBounds, {'numeric'}, ...
                {'vector', 'numel', 2, 'positive'}, 'FMAM_ODE', 'continuationOptions.quadraticStepRatioBounds');

            opts.quadraticStepRatioBounds = reshape(double(opts.quadraticStepRatioBounds),1,[]);
            if opts.quadraticStepRatioBounds(1) > opts.quadraticStepRatioBounds(2)
                error('FMAM_ODE:InvalidContinuationOptions', ...
                    'continuationOptions.quadraticStepRatioBounds must be sorted in ascending order.')
            end

            if isstring(opts.predictorMode)
                opts.predictorMode = char(opts.predictorMode);
            end
            if ~ischar(opts.predictorMode)
                error('FMAM_ODE:InvalidContinuationOptions', ...
                    'continuationOptions.predictorMode must be a string or character vector.')
            end

            opts.predictorMode = lower(opts.predictorMode);
            validModes = {'auto', 'secant', 'quadratic', 'hermite', 'constant'};
            if ~ismember(opts.predictorMode,validModes)
                error('FMAM_ODE:InvalidContinuationOptions', ...
                    'Unsupported continuationOptions.predictorMode ''%s''.', opts.predictorMode)
            end
        end

        function path = buildContinuationPath(obj)
            runtimeContext = obj.currentRuntimeContext();
            solverView = runtimeContext.solverView;
            targetCtx = runtimeContext.targetCtx;
            currentTargets = obj.initialTargetValues( ...
                solverView,runtimeContext.derived,targetCtx);
            obj.targetCurr = currentTargets;

            if obj.n_per == 0
                path = struct([]);
                return
            end

            template = struct( ...
                'startValue', 0, ...
                'rawTarget', 0, ...
                'finalValue', 0, ...
                'deltaValue', 0, ...
                'isWrapped', false, ...
                'wrapPeriod', NaN);
            path = repmat(template,1,obj.n_per);

            for i = 1:obj.n_per
                rawTarget = obj.items_perturb(i).target;
                [deltaValue,isWrapped,wrapPeriod] = fmam_target_rules('continuation_delta', ...
                    targetCtx,obj.items_perturb(i),currentTargets(i),rawTarget);
                path(i).startValue = currentTargets(i);
                path(i).rawTarget = rawTarget;
                path(i).finalValue = currentTargets(i) + deltaValue;
                path(i).deltaValue = deltaValue;
                path(i).isWrapped = isWrapped;
                path(i).wrapPeriod = wrapPeriod;
            end
        end

        function values = evaluateTargetPath(~,path,lambda)
            if isempty(path)
                values = zeros(1,0);
                return
            end

            startValues = [path.startValue];
            deltaValues = [path.deltaValue];
            values = startValues + lambda * deltaValues;
        end

        function ctx = targetRuleContext(obj,solverView,derived)
            if nargin < 2 || isempty(solverView)
                solverView = obj.solverView;
            end
            if nargin < 3 || isempty(derived)
                derived = obj.buildDerivedView(solverView);
            end
            ctx = struct( ...
                'obs', {obj.obs}, ...
                'solver', solverView, ...
                'derived', derived, ...
                'dimVar', obj.dimVar, ...
                'dimObs', obj.dimObs, ...
                'dimParams', obj.dimParams, ...
                'propertySizes', solverView.propertySizes);
        end

        function dlambdaCap = computeLambdaStepCap(obj,path)
            if isempty(path)
                dlambdaCap = 1;
                return
            end

            deltaValues = abs([path.deltaValue]);
            active = deltaValues > obj.accuracy;
            if ~any(active)
                dlambdaCap = 1;
                return
            end

            caps = obj.maxstepsize(active) ./ deltaValues(active);
            dlambdaCap = min(caps);
        end

        function dlambda = chooseInitialLambdaStep(obj,dlambdaCap)
            if ~(dlambdaCap > 0)
                dlambda = 0;
                return
            end

            if isfinite(obj.continuationOptions.initialLambdaStep) && obj.continuationOptions.initialLambdaStep > 0
                initialStep = obj.continuationOptions.initialLambdaStep;
            else
                initialStep = 1 / obj.continuationOptions.initialSteps;
            end

            dlambda = min([initialStep, dlambdaCap, obj.continuationOptions.maxLambdaStep, 1]);
        end

        function dlambda = chooseAcceptedLambdaStep(obj,acceptedStep,fitResult,dlambdaCap)
            dlambda = acceptedStep;
            if fitResult.iterations <= obj.continuationOptions.easyIterations
                dlambda = dlambda * obj.continuationOptions.growthFactor;
            elseif fitResult.iterations >= obj.continuationOptions.slowIterations
                dlambda = dlambda * obj.continuationOptions.shrinkFactor;
            end

            dlambda = obj.clampLambdaStep(dlambda,dlambdaCap,1);
        end

        function dlambda = chooseRejectedLambdaStep(obj,rejectedStep,dlambdaCap)
            dlambda = obj.clampLambdaStep( ...
                rejectedStep * obj.continuationOptions.backtrackingFactor, dlambdaCap, 1);
        end

        function dlambda = clampLambdaStep(obj,requested,dlambdaCap,remaining)
            values = [requested, dlambdaCap, obj.continuationOptions.maxLambdaStep, remaining];
            values = values(isfinite(values));
            if isempty(values)
                dlambda = 0;
            else
                dlambda = max(0,min(values));
            end
        end

        function acceptedState = captureAcceptedContinuationState(obj,lambda,dlambdaAccepted,fitResult,predictorInfo)
            if nargin < 5 || isempty(predictorInfo)
                predictorInfo = struct();
            end
            if ~isfield(predictorInfo,'used') || isempty(predictorInfo.used)
                predictorInfo.used = 'constant';
            end
            if ~isfield(predictorInfo,'fallbackCount') || isempty(predictorInfo.fallbackCount)
                predictorInfo.fallbackCount = 0;
            end
            if ~isfield(predictorInfo,'message') || isempty(predictorInfo.message)
                predictorInfo.message = '';
            end

            snapshot = obj.snapshotSolverView();
            [A,~,meta] = obj.linearizeNewtonProblem();
            diagnostics = struct( ...
                'iterations', NaN, ...
                'message', '', ...
                'usedLM', false, ...
                'finalConditionEstimate', NaN, ...
                'directConditionEstimate', NaN);

            if nargin >= 4 && ~isempty(fitResult) && isstruct(fitResult)
                if isfield(fitResult,'iterations')
                    diagnostics.iterations = fitResult.iterations;
                end
                if isfield(fitResult,'message')
                    diagnostics.message = fitResult.message;
                end
                if isfield(fitResult,'history') && ~isempty(fitResult.history)
                    diagnostics.usedLM = any(strcmp({fitResult.history.solver},'lm'));
                    diagnostics.finalConditionEstimate = fitResult.history(end).conditionEstimate;
                    if isfield(fitResult.history,'directConditionEstimate')
                        diagnostics.directConditionEstimate = ...
                            fitResult.history(end).directConditionEstimate;
                    end
                end
            end

            acceptedState = struct( ...
                'point', struct( ...
                    'lambda', lambda, ...
                    'packedUnknown', obj.packUnknownState(meta,snapshot), ...
                    'dlambdaAccepted', dlambdaAccepted, ...
                    'newtonDiagnostics', diagnostics, ...
                    'tangent', obj.computeContinuationTangent(A), ...
                    'predictorUsed', predictorInfo.used, ...
                    'predictorFallbackCount', predictorInfo.fallbackCount, ...
                    'predictorMessage', predictorInfo.message), ...
                'snapshot', snapshot, ...
                'meta', meta, ...
                'targetCurr', obj.targetCurr, ...
                'stepsize', obj.stepsize);
        end

        function restoreAcceptedContinuationState(obj,acceptedState)
            obj.restoreSolverView(acceptedState.snapshot);
            obj.targetCurr = acceptedState.targetCurr;
            if isfield(acceptedState,'stepsize')
                obj.stepsize = acceptedState.stepsize;
            end
        end

        function activateTrialContinuationState(obj,acceptedState,trialTarget,trialSnapshot)
            try
                obj.targetCurr = trialTarget;
                obj.stepsize = trialTarget - acceptedState.targetCurr;
                obj.restoreSolverView(trialSnapshot);
            catch ME
                obj.restoreAcceptedContinuationState(acceptedState);
                rethrow(ME)
            end
        end

        function packed = packUnknownState(~,meta,snapshot)
            totalSize = 0;
            for i = 1:numel(meta.unknowns)
                totalSize = max(totalSize,max(meta.indexMap.(meta.unknowns{i})(:)));
            end

            packed = zeros(totalSize,1);
            for i = 1:numel(meta.unknowns)
                propname = meta.unknowns{i};
                packed(meta.indexMap.(propname)(:)) = snapshot.(propname)(:);
            end
        end

        function issue = convergedStatePsiIssue(obj)
            if ~obj.checkPsiNonnegative
                issue = struct('isValid', true, 'identifier', '', 'message', '');
                return
            end
            issue = fmam_state_ops.positivePsiIssue(obj.currentPsiValues());
        end

        function [failureDiagnostics, failCount, baseFailureCount, dlambda, stopNow] = ...
                rejectContinuationTrial(obj,failureDiagnostics,failCount,baseFailureCount, ...
                acceptedState,count,lambda,lambdaTrial, ...
                dlambdaTrial,dlambdaCap,fitResult,rejectionMessage)
            obj.restoreAcceptedContinuationState(acceptedState);

            failureEntry = struct( ...
                'trialIteration', count, ...
                'lambdaTrial', lambdaTrial, ...
                'dlambdaTrial', dlambdaTrial, ...
                'fitResult', fitResult);
            failureDiagnostics(end + 1) = failureEntry;
            failCount = failCount + 1;
            baseFailureCount = baseFailureCount + 1;
            dlambda = obj.chooseRejectedLambdaStep(dlambdaTrial,dlambdaCap);
            stopNow = failCount > obj.continuationOptions.maxFailures;

            obj.logMessage('REJECT', ...
                ['iter=%d rejected at lambdaTrial=%.6f: message="%s", fail=%d/%d, ' ...
                'next dlambda=%.6e'], ...
                count, lambdaTrial, rejectionMessage, failCount, ...
                obj.continuationOptions.maxFailures, dlambda);

            if stopNow
                obj.logMessage('STOP', ...
                    ['failCount=%d exceeded maxFailures=%d at lambda=%.6f; ' ...
                    'terminate continuation.'], ...
                    failCount, obj.continuationOptions.maxFailures, lambda);
            end
        end

        function tangent = computeContinuationTangent(obj,A)
            tangent = struct( ...
                'available', false, ...
                'vector', [], ...
                'solver', '', ...
                'conditionEstimate', NaN, ...
                'lambda', NaN, ...
                'message', '');

            if isempty(obj.targetPath)
                tangent.message = 'continuation path is not initialized';
                return
            end

            rhs = zeros(size(A,1),1);
            rowIdx = obj.targetEquationRowIndices(size(A,1));
            if isempty(rowIdx)
                tangent.available = true;
                tangent.vector = zeros(size(A,2),1);
                tangent.message = 'no active continuation targets';
                return
            end
            rhs(rowIdx) = reshape([obj.targetPath.deltaValue],[],1);

            solveResult = solve_regularized_linear_system( ...
                A,rhs,obj.linearSolverOptions(struct()),'best_effort');
            tangent.available = solveResult.success;
            tangent.vector = solveResult.solution;
            tangent.solver = solveResult.solver;
            tangent.conditionEstimate = solveResult.conditionEstimate;
            tangent.lambda = solveResult.lambda;
            tangent.message = solveResult.message;
        end

        function rowIdx = targetEquationRowIndices(obj,numRows)
            if obj.n_per == 0
                rowIdx = zeros(1,0);
                return
            end

            if obj.isPsiUpdated
                trailingConstraints = max(0,2 * obj.truncationOrder - 1);
            else
                trailingConstraints = max(0,2 * obj.truncationOrder - 2);
            end

            rowStart = numRows - trailingConstraints - obj.n_per + 1;
            rowIdx = rowStart:(rowStart + obj.n_per - 1);
        end

        function info = preparePredictorTrial(obj,acceptedState,acceptedHistory,lambdaTrial,dlambdaTrial,baseFailureCount)
            info = struct( ...
                'used', 'constant', ...
                'fallbackCount', 0, ...
                'message', 'constant predictor', ...
                'snapshot', acceptedState.snapshot);

            predictorOrder = obj.predictorOrder(acceptedHistory,dlambdaTrial,baseFailureCount);
            fallbackCount = 0;
            lastMessage = 'constant predictor';
            for i = 1:numel(predictorOrder)
                predictorMode = predictorOrder{i};
                candidate = obj.computePredictorCandidate(acceptedHistory,lambdaTrial,predictorMode);
                if ~candidate.available
                    lastMessage = candidate.message;
                    fallbackCount = fallbackCount + 1;
                    continue
                end

                validation = obj.validatePredictorCandidate( ...
                    acceptedState,acceptedHistory,candidate.packedUnknown,predictorMode);
                if validation.valid
                    info.used = predictorMode;
                    info.fallbackCount = fallbackCount;
                    info.message = validation.message;
                    info.snapshot = validation.snapshot;
                    return
                end

                lastMessage = validation.message;
                fallbackCount = fallbackCount + 1;
            end

            info.message = lastMessage;
        end

        function order = predictorOrder(obj,acceptedHistory,dlambdaTrial,baseFailureCount)
            if baseFailureCount >= 2
                order = {'constant'};
                return
            end

            predictorMode = obj.continuationOptions.predictorMode;
            if strcmp(predictorMode,'auto')
                predictorMode = obj.selectAutoPredictor(acceptedHistory,dlambdaTrial);
            end

            order = FMAM_ODE.predictorFallbackOrder(predictorMode);
            if baseFailureCount >= 1
                order(strcmp(order,'quadratic')) = [];
                if isempty(order)
                    order = {'constant'};
                elseif ~strcmp(order{end},'constant')
                    order{end + 1} = 'constant';
                end
            end
        end

        function predictorMode = selectAutoPredictor(obj,acceptedHistory,dlambdaTrial)
            predictorMode = 'constant';
            numAccepted = numel(acceptedHistory);
            if numAccepted >= 3
                older = acceptedHistory(end - 2);
                prev = acceptedHistory(end - 1);
                last = acceptedHistory(end);
                d1 = prev.lambda - older.lambda;
                d2 = last.lambda - prev.lambda;

                if d1 > 0 && d2 > 0
                    ratio = d2 / d1;
                    slope1 = (prev.packedUnknown - older.packedUnknown) / d1;
                    slope2 = (last.packedUnknown - prev.packedUnknown) / d2;
                    curvature = norm(slope2 - slope1,inf) / max(norm(slope2,inf),eps);
                    bounds = obj.continuationOptions.quadraticStepRatioBounds;
                    if ratio >= bounds(1) && ratio <= bounds(2) && ...
                            curvature <= obj.continuationOptions.quadraticCurvatureThreshold
                        predictorMode = 'quadratic';
                        return
                    end
                end
            end

            if numAccepted >= 2
                prev = acceptedHistory(end - 1);
                last = acceptedHistory(end);
                lastAcceptedStep = last.lambda - prev.lambda;
                if lastAcceptedStep > 0 && prev.tangent.available && last.tangent.available && ...
                        dlambdaTrial / lastAcceptedStep <= obj.continuationOptions.hermiteMaxExtrapolationRatio
                    predictorMode = 'hermite';
                else
                    predictorMode = 'secant';
                end
            end
        end

        function candidate = computePredictorCandidate(~,acceptedHistory,lambdaTrial,predictorMode)
            basePoint = acceptedHistory(end);
            candidate = struct('available', false, 'packedUnknown', basePoint.packedUnknown, 'message', '');

            switch predictorMode
                case 'constant'
                    candidate.available = true;
                    candidate.message = 'constant predictor';

                case 'secant'
                    if numel(acceptedHistory) < 2
                        candidate.message = 'secant predictor needs at least two accepted points';
                        return
                    end
                    prev = acceptedHistory(end - 1);
                    last = acceptedHistory(end);
                    denom = last.lambda - prev.lambda;
                    if abs(denom) <= eps
                        candidate.message = 'secant predictor requires distinct continuation parameters';
                        return
                    end
                    factor = (lambdaTrial - last.lambda) / denom;
                    candidate.packedUnknown = last.packedUnknown + factor * (last.packedUnknown - prev.packedUnknown);
                    candidate.available = true;
                    candidate.message = 'secant predictor';

                case 'quadratic'
                    if numel(acceptedHistory) < 3
                        candidate.message = 'quadratic predictor needs at least three accepted points';
                        return
                    end
                    p0 = acceptedHistory(end - 2);
                    p1 = acceptedHistory(end - 1);
                    p2 = acceptedHistory(end);
                    lambdas = [p0.lambda, p1.lambda, p2.lambda];
                    if min(abs(diff(lambdas))) <= eps || abs(lambdas(3) - lambdas(1)) <= eps
                        candidate.message = 'quadratic predictor requires distinct continuation parameters';
                        return
                    end
                    L0 = (lambdaTrial - lambdas(2)) * (lambdaTrial - lambdas(3)) / ...
                        ((lambdas(1) - lambdas(2)) * (lambdas(1) - lambdas(3)));
                    L1 = (lambdaTrial - lambdas(1)) * (lambdaTrial - lambdas(3)) / ...
                        ((lambdas(2) - lambdas(1)) * (lambdas(2) - lambdas(3)));
                    L2 = (lambdaTrial - lambdas(1)) * (lambdaTrial - lambdas(2)) / ...
                        ((lambdas(3) - lambdas(1)) * (lambdas(3) - lambdas(2)));
                    candidate.packedUnknown = L0 * p0.packedUnknown + L1 * p1.packedUnknown + L2 * p2.packedUnknown;
                    candidate.available = true;
                    candidate.message = 'quadratic predictor';

                case 'hermite'
                    if numel(acceptedHistory) < 2
                        candidate.message = 'Hermite predictor needs at least two accepted points';
                        return
                    end
                    p0 = acceptedHistory(end - 1);
                    p1 = acceptedHistory(end);
                    if ~p0.tangent.available || ~p1.tangent.available
                        candidate.message = 'Hermite predictor requires finite endpoint tangents';
                        return
                    end
                    h = p1.lambda - p0.lambda;
                    if h <= eps
                        candidate.message = 'Hermite predictor requires distinct continuation parameters';
                        return
                    end
                    t = (lambdaTrial - p0.lambda) / h;
                    h00 = 2 * t^3 - 3 * t^2 + 1;
                    h10 = t^3 - 2 * t^2 + t;
                    h01 = -2 * t^3 + 3 * t^2;
                    h11 = t^3 - t^2;
                    candidate.packedUnknown = ...
                        h00 * p0.packedUnknown + ...
                        h10 * h * p0.tangent.vector + ...
                        h01 * p1.packedUnknown + ...
                        h11 * h * p1.tangent.vector;
                    candidate.available = true;
                    candidate.message = 'Hermite predictor';

                otherwise
                    candidate.message = sprintf('Unsupported predictor ''%s''.', predictorMode);
            end
        end

        function validation = validatePredictorCandidate(obj,acceptedState,acceptedHistory,packedUnknown,predictorMode)
            basePoint = acceptedState.point;
            validation = struct('valid', false, 'snapshot', acceptedState.snapshot, 'message', '');

            if numel(packedUnknown) ~= numel(basePoint.packedUnknown) || any(~isfinite(packedUnknown))
                validation.message = sprintf('%s predictor produced a non-finite iterate.', predictorMode);
                return
            end

            if numel(acceptedHistory) >= 2
                referenceStep = norm(acceptedHistory(end).packedUnknown - acceptedHistory(end - 1).packedUnknown,inf);
            else
                referenceStep = 0;
            end
            predictorJump = norm(packedUnknown - basePoint.packedUnknown,inf);
            jumpLimit = obj.continuationOptions.predictorStepGrowthLimit * max(1,referenceStep);
            if predictorJump > jumpLimit
                validation.message = sprintf('%s predictor exceeded the step-growth limit.', predictorMode);
                return
            end

            obj.restoreSolverView(acceptedState.snapshot);
            try
                obj.applyUnknownIncrement(acceptedState.meta,packedUnknown - basePoint.packedUnknown,1);
                predictedSnapshot = obj.snapshotSolverView();
            catch ME
                obj.restoreSolverView(acceptedState.snapshot);
                validation.message = sprintf('%s predictor was rejected: %s', ...
                    predictorMode, FMAM_ODE.exceptionSummary(ME));
                return
            end

            obj.restoreSolverView(acceptedState.snapshot);
            validation.valid = true;
            validation.snapshot = predictedSnapshot;
            validation.message = sprintf('%s predictor accepted', predictorMode);
        end

        function synchronizePsiReferenceMode(obj)
            if ~obj.isPsiUpdated && (isempty(obj.p_Psi_init) || isempty(obj.q_Psi_init))
                obj.captureFrozenPsiReference();
            end
        end

        function applyPsiUpdateMode(obj,value)
            validateattributes(value, {'numeric', 'logical'}, ...
                {'scalar'}, 'FMAM_ODE', 'isPsiUpdated');
            nextMode = logical(value);
            previousMode = obj.psiUpdateMode;
            obj.psiUpdateMode = nextMode;

            if ~nextMode
                if previousMode || isempty(obj.p_Psi_init) || isempty(obj.q_Psi_init)
                    obj.captureFrozenPsiReference();
                end
            end
        end

        function captureFrozenPsiReference(obj)
            obj.p_Psi_init = obj.solverView.p_Psi;
            obj.q_Psi_init = obj.solverView.q_Psi;
        end

        function validateStateDimensions(obj,solverView)
            if ~iscell(obj.sys)
                error('system must be a cell array of ODE right-hand-side functions.')
            end
            if ~isempty(obj.obs) && ~iscell(obj.obs)
                error('observables must be provided as a cell array.')
            end
            if obj.dimVar ~= solverView.dimSys
                error('The dimension of system does not match solverView.dimSys.')
            end
            if obj.dimObs ~= solverView.dimObs
                error('The number of observables does not match solverView.dimObs.')
            end
        end

        function validateDerivativeCache(obj)
            if isempty(obj.derivatives)
                error('FMAM_ODE:MissingDerivatives', ...
                    ['FMAM_ODE requires precomputed derivatives. Pass a derivatives struct ' ...
                    'during construction, for example ''derivatives'', build_symbolic_derivatives(...).'])
            end
            if ~isstruct(obj.derivatives) || ~all(isfield(obj.derivatives,{'var','obs','obs2'}))
                error('FMAM_ODE:InvalidDerivatives', ...
                    'derivatives must be a struct with fields var, obs, and obs2.')
            end

            obj.validateDerivativeEntries(obj.derivatives.var,[obj.dimVar,obj.dimVar + obj.dimParams],'derivatives.var')

            if obj.dimObs == 0
                if ~isempty(obj.derivatives.obs) || ~isempty(obj.derivatives.obs2)
                    error('FMAM_ODE:InvalidDerivatives', ...
                        'derivatives.obs and derivatives.obs2 must be empty when no observables are provided.')
                end
                return
            end

            obj.validateDerivativeEntries(obj.derivatives.obs,[obj.dimObs,obj.dimVar],'derivatives.obs')
            obj.validateDerivativeEntries(obj.derivatives.obs2,[obj.dimObs,obj.dimVar,obj.dimVar],'derivatives.obs2')
        end

        function validateDerivativeEntries(~,entries,expectedSize,entryName)
            if ~isstruct(entries) || ~isequal(size(entries),expectedSize) || ~all(isfield(entries,'function'))
                error('FMAM_ODE:InvalidDerivatives', ...
                    '%s must be a struct array of size %s with a field named function.', ...
                    entryName, mat2str(expectedSize))
            end

            handles = {entries.function};
            if ~all(cellfun(@(f) isa(f,'function_handle'), handles))
                error('FMAM_ODE:InvalidDerivatives', ...
                    'Every entry in %s must store a function handle in its function field.', entryName)
            end
        end

        function appendAcceptedLog(obj,acceptedPoint,fitResult)
            runtimeContext = obj.currentRuntimeContext();
            currentView = runtimeContext.solverView;
            derived = runtimeContext.derived;
            log_curr = struct('params',currentView.params);
            targetCtx = runtimeContext.targetCtx;

            for i = 1:obj.n_per
                item_per = obj.items_perturb(i);
                propname = item_per.prop;
                idxLabel = sprintf('%g_',item_per.idx(:));
                idxLabel = regexprep(idxLabel(1:end-1),'[^A-Za-z0-9_]','_');
                logName = [propname '_idx_' idxLabel];
                log_curr.(logName) = obj.targetCurr(i);
                log_curr.(['accepted_' logName]) = obj.currentTargetValue( ...
                    item_per,currentView,derived,targetCtx);
            end

            if nargin >= 2 && isstruct(acceptedPoint)
                log_curr.lambda = acceptedPoint.lambda;
                log_curr.dlambdaAccepted = acceptedPoint.dlambdaAccepted;
                log_curr.predictorUsed = acceptedPoint.predictorUsed;
                log_curr.predictorFallbackCount = acceptedPoint.predictorFallbackCount;
            else
                log_curr.lambda = NaN;
                log_curr.dlambdaAccepted = NaN;
                log_curr.predictorUsed = '';
                log_curr.predictorFallbackCount = 0;
            end

            if nargin >= 3 && isstruct(fitResult)
                log_curr.newtonIterations = fitResult.iterations;
                log_curr.message = fitResult.message;
                if isfield(fitResult,'history') && ~isempty(fitResult.history)
                    log_curr.finalConditionEstimate = fitResult.history(end).conditionEstimate;
                    if isfield(fitResult.history,'directConditionEstimate')
                        log_curr.directConditionEstimate = ...
                            fitResult.history(end).directConditionEstimate;
                    else
                        log_curr.directConditionEstimate = NaN;
                    end
                    log_curr.usedLM = any(strcmp({fitResult.history.solver},'lm'));
                else
                    log_curr.finalConditionEstimate = NaN;
                    log_curr.directConditionEstimate = NaN;
                    log_curr.usedLM = false;
                end
            else
                log_curr.newtonIterations = NaN;
                log_curr.message = '';
                log_curr.finalConditionEstimate = NaN;
                log_curr.directConditionEstimate = NaN;
                log_curr.usedLM = false;
            end

            stateMeasureIdx = min(2, obj.dimVar);
            log_curr.period = derived.period;
            log_curr.amplitude = derived.varAmp(stateMeasureIdx);
            log_curr.yMax = derived.varMax(stateMeasureIdx);
            log_curr.yMin = derived.varMin(stateMeasureIdx);

            if isempty(obj.logs)
                obj.logs = log_curr;
            else
                obj.logs = [obj.logs log_curr];
            end
        end

        function [itemsPer,controlled,maxstepsize,accuracy] = validateModulationSetup(obj,itemsPer,controlled,maxstepsize,accuracy)
            if ~isstruct(itemsPer)
                error('items_perturb must be a struct array.')
            end
            itemsPer = reshape(itemsPer,1,[]);
            controlled = reshape(controlled,1,[]);
            nPer = numel(itemsPer);

            if ~isempty(itemsPer) && ~all(isfield(itemsPer,{'prop','idx','target'}))
                error('Each items_perturb entry must contain the fields ''prop'', ''idx'', and ''target''.')
            end
            if any(controlled ~= floor(controlled)) || any(controlled < 1) || any(controlled > obj.dimParams)
                error('items_controlled must contain unique parameter indices within the parameter range.')
            end
            if numel(unique(controlled)) ~= numel(controlled)
                error('items_controlled must not contain duplicate parameter indices.')
            end
            if numel(controlled) ~= nPer
                error('numel(items_controlled) must equal numel(items_perturb).')
            end

            maxstepsize = obj.normalizeContinuationArgument(maxstepsize,nPer,'maxstepsize');
            accuracy = obj.normalizeContinuationArgument(accuracy,nPer,'accuracy');

            for i = 1:nPer
                itemsPer(i).prop = obj.canonicalTargetProperty(itemsPer(i).prop);
                if ~isnumeric(itemsPer(i).target) || ~isscalar(itemsPer(i).target) || ~isfinite(itemsPer(i).target)
                    error('Each modulation target must be a finite numeric scalar.')
                end
                obj.validateTargetIndex(itemsPer(i));
                if strcmp(itemsPer(i).prop,'params') && ~ismember(itemsPer(i).idx,controlled)
                    error('Direct parameter targets must also appear in items_controlled.')
                end
            end
        end

        function values = normalizeContinuationArgument(~,values,nPer,argName)
            if nPer == 0
                values = zeros(1,0);
                return
            end
            if isempty(values)
                if strcmp(argName,'maxstepsize')
                    values = inf(1,nPer);
                    return
                end
                error('%s must not be empty.',argName)
            end
            if isscalar(values)
                values = repmat(double(values),1,nPer);
            else
                values = reshape(double(values),1,[]);
            end
            if numel(values) ~= nPer
                error('%s must be a scalar or a vector with one entry per modulation target.',argName)
            end
            if any(~isfinite(values)) || any(values < 0)
                error('%s must contain finite nonnegative values.',argName)
            end
        end

        function propname = canonicalTargetProperty(~,propname)
            propname = fmam_target_canonicalize(propname);
        end

        function validateTargetIndex(obj,item)
            fmam_target_rules('validate_item',obj.targetRuleContext(),item);
        end

        function targets = initialTargetValues(obj,solverView,derived,targetCtx)
            if nargin < 2
                solverView = [];
            end
            if nargin < 3
                derived = [];
            end
            if nargin < 4 || isempty(targetCtx)
                runtimeContext = obj.currentRuntimeContext(solverView,derived);
                solverView = runtimeContext.solverView;
                derived = runtimeContext.derived;
                targetCtx = runtimeContext.targetCtx;
            end
            targets = zeros(1,obj.n_per);

            for i = 1:obj.n_per
                targets(i) = obj.currentTargetValue(obj.items_perturb(i),solverView,derived,targetCtx);
            end
        end

        function [needMax,needMin] = targetNeedsVariableExtrema(obj,varIdx,targetCtx)
            if nargin < 3 || isempty(targetCtx)
                targetCtx = obj.currentRuntimeContext().targetCtx;
            end
            [needMax,needMin] = fmam_target_rules('needs_variable_extrema', ...
                targetCtx,obj.items_perturb,varIdx);
        end

        function [needMax,needMin] = targetNeedsObservableExtrema(obj,obsIdx,targetCtx)
            if nargin < 3 || isempty(targetCtx)
                targetCtx = obj.currentRuntimeContext().targetCtx;
            end
            [needMax,needMin] = fmam_target_rules('needs_observable_extrema', ...
                targetCtx,obj.items_perturb,obsIdx);
        end

        function value = currentTargetValue(obj,item,solverView,derived,targetCtx)
            if nargin < 3
                solverView = [];
            end
            if nargin < 4
                derived = [];
            end
            if nargin < 5 || isempty(targetCtx)
                runtimeContext = obj.currentRuntimeContext(solverView,derived);
                solverView = runtimeContext.solverView;
                targetCtx = runtimeContext.targetCtx;
            end
            value = fmam_target_rules('current_value',targetCtx, ...
                item,solverView);
        end

        function logNewtonSummary(obj,fitResult)
            if ~obj.verbose
                return
            end

            obj.logMessage('NEWTON', ...
                ['summary: converged=%d, iterations=%d, maxErr=%s, residual2=%s, ' ...
                'message="%s"'], ...
                fitResult.converged, fitResult.iterations, ...
                obj.formatScalar(fitResult.scalarError), ...
                obj.formatScalar(fitResult.linearResidualNorm), fitResult.message);

            history = fitResult.history;
            for i = 1:numel(history)
                h = history(i);
                obj.logMessage('NEWTON', ...
                    ['iter=%d, solver=%s, accepted=%d, stepInf=%s, residual2=%s, ' ...
                    'maxErr=%s, cond=%s, directCond=%s, lambda=%s'], ...
                    h.iteration, h.solver, h.accepted, ...
                    obj.formatScalar(h.stepNorm), obj.formatScalar(h.residualNorm), ...
                    obj.formatScalar(h.maxError), ...
                    obj.formatScalar(h.conditionEstimate), ...
                    obj.formatScalar(h.directConditionEstimate), ...
                    obj.formatScalar(h.lambda));
            end
        end

        function tf = shouldStopForCondition(obj,directCondition)
            tf = logical(obj.continuationOptions.conditionStopEnabled) && ...
                isfinite(directCondition) && ...
                directCondition < obj.continuationOptions.conditionStopRcond;
        end

        function value = extractDirectConditionEstimate(~,fitResult)
            value = NaN;
            if ~isstruct(fitResult) || ~isfield(fitResult,'history') || isempty(fitResult.history)
                return
            end
            if isfield(fitResult.history,'directConditionEstimate')
                value = fitResult.history(end).directConditionEstimate;
            end
        end

        function logMessage(obj,stage,fmt,varargin)
            if ~obj.verbose
                return
            end
            message = sprintf(fmt,varargin{:});
            fprintf('[FMAM_ODE][%s] %s\n',stage,message);
        end

        function text = formatScalar(~,value)
            if isempty(value) || ~isscalar(value) || ~isnumeric(value) || ~isfinite(value)
                text = 'NaN';
            else
                text = sprintf('%.6e',value);
            end
        end

        function result = translateNewtonResult(~,solveResult,errorFieldName)
            result = struct();
            result.converged = solveResult.converged;
            result.iterations = solveResult.iterations;
            result.(errorFieldName) = solveResult.errorVec;
            result.message = solveResult.message;
            result.history = solveResult.history;
            result.linearResidualNorm = solveResult.linearResidualNorm;
            result.linearResidual = solveResult.linearResidual;
            result.stepAccepted = solveResult.stepAccepted;
            result.scalarError = solveResult.scalarError;
            result.objective = solveResult.objective;
        end
    end

    methods (Static, Access = private)
        function base = mergeOptionStruct(base,overrides,argName,errorId)
            if nargin < 3
                argName = 'options';
            end
            if nargin < 4 || isempty(errorId)
                errorId = 'FMAM_ODE:UnknownOption';
            end
            if nargin < 2 || isempty(overrides)
                return
            end
            if ~isstruct(overrides)
                error('FMAM_ODE:InvalidNewtonOptions','%s must be a struct.',argName)
            end

            unknown = setdiff(fieldnames(overrides), fieldnames(base));
            if ~isempty(unknown)
                error(errorId, 'Unknown %s field(s): %s.', argName, strjoin(unknown, ', '));
            end

            names = fieldnames(overrides);
            for i = 1:numel(names)
                base.(names{i}) = overrides.(names{i});
            end
        end

        function opts = defaultContinuationOptions()
            opts = struct();
            opts.initialSteps = 8;
            opts.initialLambdaStep = NaN;
            opts.minLambdaStep = 1e-4;
            opts.maxLambdaStep = 1;
            opts.growthFactor = 1.1;
            opts.shrinkFactor = 0.75;
            opts.backtrackingFactor = 0.5;
            opts.easyIterations = 4;
            opts.slowIterations = 10;
            opts.maxFailures = 8;
            opts.conditionStopEnabled = false;
            opts.conditionStopRcond = 1e-8;
            opts.predictorMode = 'auto';
            opts.quadraticCurvatureThreshold = 0.25;
            opts.quadraticStepRatioBounds = [0.5, 2.0];
            opts.hermiteMaxExtrapolationRatio = 1.5;
            opts.predictorStepGrowthLimit = 3.0;
        end

        function order = predictorFallbackOrder(mode)
            switch mode
                case 'quadratic'
                    order = {'quadratic', 'hermite', 'secant', 'constant'};
                case 'hermite'
                    order = {'hermite', 'secant', 'constant'};
                case 'secant'
                    order = {'secant', 'constant'};
                case 'constant'
                    order = {'constant'};
                otherwise
                    error('FMAM_ODE:InvalidPredictorMode', ...
                        'Unsupported predictor mode ''%s''.', mode)
            end
        end

        function delta = wrapPeriodicDifference(delta,period)
            if ~(isfinite(period) && period > 0)
                return
            end

            delta = mod(delta + 0.5 * period, period) - 0.5 * period;
        end

        function summary = exceptionSummary(ME)
            if ~isempty(ME.identifier)
                summary = ME.identifier;
            else
                summary = ME.message;
            end
        end

        function view = applyIncrementToSolverView(view,meta,increments,scale)
            if nargin < 4 || isempty(scale)
                scale = 1;
            end
            if ~isstruct(meta) || ~isfield(meta,'unknowns') || ~isfield(meta,'indexMap')
                error('FMAM_ODE:InvalidSolverMeta', ...
                    'meta must contain unknowns and indexMap fields.')
            end

            scaledIncrement = scale * increments;
            for i = 1:numel(meta.unknowns)
                propname = meta.unknowns{i};
                if ~isfield(meta.indexMap,propname)
                    error('FMAM_ODE:InvalidSolverMeta', ...
                        'indexMap is missing the unknown block ''%s''.', propname)
                end
                idxProp = meta.indexMap.(propname);
                blockValue = view.(propname);
                columns = size(idxProp,2);
                blockIncrement = reshape(scaledIncrement(idxProp(:)),[],columns);
                if ~isequal(size(blockValue),size(blockIncrement))
                    blockIncrement = reshape(blockIncrement,size(blockValue));
                end
                view.(propname) = blockValue + blockIncrement;
            end
        end

    end

    methods (Access = private, Static)
        function opts = parseConstructorOptions(varargin)
            if mod(numel(varargin),2) ~= 0
                error('FMAM_ODE:InvalidConstructorOption', ...
                    'Constructor options must be supplied as name-value pairs.')
            end

            opts = struct( ...
                'derivatives', [], ...
                'newtonOptions', struct(), ...
                'continuationOptions', struct(), ...
                'isPsiUpdated', false, ...
                'checkPsiNonnegative', true, ...
                'needLog', false, ...
                'verbose', true, ...
                'errBound', 1e-8, ...
                'Lconst', 500);

            for k = 1:2:numel(varargin)
                name = varargin{k};
                value = varargin{k+1};
                if isstring(name) && isscalar(name)
                    name = char(name);
                end
                if ~ischar(name)
                    error('FMAM_ODE:InvalidConstructorOption', ...
                        'Constructor option names must be character vectors or strings.')
                end

                switch lower(name)
                    case 'derivatives'
                        opts.derivatives = value;
                    case 'newtonoptions'
                        if ~isstruct(value)
                            error('FMAM_ODE:InvalidConstructorOption', ...
                                'Constructor option ''newtonOptions'' must be a struct.')
                        end
                        opts.newtonOptions = value;
                    case 'continuationoptions'
                        if ~isstruct(value)
                            error('FMAM_ODE:InvalidConstructorOption', ...
                                'Constructor option ''continuationOptions'' must be a struct.')
                        end
                        opts.continuationOptions = value;
                    case 'ispsiupdated'
                        validateattributes(value, {'numeric', 'logical'}, ...
                            {'scalar'}, 'FMAM_ODE', 'isPsiUpdated');
                        opts.isPsiUpdated = logical(value);
                    case 'checkpsinonnegative'
                        validateattributes(value, {'numeric', 'logical'}, ...
                            {'scalar'}, 'FMAM_ODE', 'checkPsiNonnegative');
                        opts.checkPsiNonnegative = logical(value);
                    case 'needlog'
                        validateattributes(value, {'numeric', 'logical'}, ...
                            {'scalar'}, 'FMAM_ODE', 'needLog');
                        opts.needLog = logical(value);
                    case 'verbose'
                        validateattributes(value, {'numeric', 'logical'}, ...
                            {'scalar'}, 'FMAM_ODE', 'verbose');
                        opts.verbose = logical(value);
                    case 'errbound'
                        validateattributes(value, {'numeric'}, ...
                            {'scalar', 'positive', 'finite'}, 'FMAM_ODE', 'errBound');
                        opts.errBound = double(value);
                    case 'lconst'
                        validateattributes(value, {'numeric'}, ...
                            {'scalar', 'integer', 'positive', 'finite'}, 'FMAM_ODE', 'Lconst');
                        opts.Lconst = double(value);
                    otherwise
                        error('FMAM_ODE:InvalidConstructorOption', ...
                            'Unsupported constructor option ''%s''.', name)
                end
            end
        end
    end

    methods(Static)
        function [vc,vs] = Vec_CS(phi, M, L)
        % Vec_CS  Build [cos(0*phi)...cos(M*phi)] and [sin(1*phi)...sin(M*phi)]
        % phi : L-by-1 (or 1-by-L) vector of angles (radians)
        % M   : highest harmonic
        % L   : (optional) asserted length; ignored if empty
        %
        % Output:
        %   vc : L-by-(M+1), columns k=0..M  -> cos(k*phi)
        %   vs : L-by-M,     columns k=1..M  -> sin(k*phi)
        
            phi = phi(:);                 % ensure column
            if nargin>=3 && ~isempty(L)
                assert(numel(phi)==L, 'L must match length(phi).');
            end
            L = numel(phi);
        
            kC = 0:M;                      % 1-by-(M+1)
            kS = 1:M;                      % 1-by-M
        
            % Use outer products via implicit expansion (no big loops in MATLAB)
            vc = cos(phi * kC);            % L x (M+1)
            if M>0
                vs = sin(phi * kS);        % L x M
            else
                vs = zeros(L,0,'like',phi);
            end
        end


        function output = Trintegration(p,q,phi_L,phi_U)
            % Compute the integration of a trigonometric polynomial with coefficients
            % saved in p and q. phi_L and phi_U are the lower bound and upper bound of
            % the integral, respectively.
            M = size(q,1);
            [vc_U,vs_U] = FMAM_ODE.Vec_CS(phi_U,M,1);
            [vc_L,vs_L] = FMAM_ODE.Vec_CS(phi_L,M,1);
            U = p(1)*phi_U+vs_U.*((1:M).^(-1))*p(2:end)-vc_U(:,2:end).*((1:M).^(-1))*q;
            L = p(1)*phi_L+vs_L.*((1:M).^(-1))*p(2:end)-vc_L(:,2:end).*((1:M).^(-1))*q;
            output = U-L;
        end

        function [coe_para,coe_p_Psi,coe_q_Psi,coe_p_var,coe_q_var] = ...
         delta_coe_system(System,Derivative,tphi,TS_variable,parameters,...
         vcPsi,vsPsi,vcVar,vsVar,idx_eq)
            MVar = size(vsVar,2);
            L = size(TS_variable,1);
            N = size(TS_variable,2);
            m = size(parameters,2);

            coe_para = zeros(L,m);

            F = System{1,idx_eq};
        
            for i = 1:m
                coe_para(:,i) = -tphi.*Derivative(idx_eq,N+i).function(TS_variable,parameters);
            end
            coe_p_Psi = -F(TS_variable,parameters).*vcPsi;
            coe_q_Psi = -F(TS_variable,parameters).*vsPsi;

            A = zeros(L,N);
            for i = 1:N 
                A(:,i) = -tphi.*Derivative(idx_eq,i).function(TS_variable,parameters);
            end

            coe_p_var = zeros(L,N*(MVar+1));
            coe_q_var = zeros(L,N*MVar);

            DVS = [zeros(size(vsVar,1),1),vsVar.*(1:MVar)];
            DVC = vcVar(:,2:end).*(1:MVar);
            for i = 1:N
                idx_p = (i-1)*(MVar+1) + (1:MVar+1);
                idx_q = (i-1)*MVar + (1:MVar);
                coe_p_var(:,idx_p) = A(:,i).*vcVar;
                coe_q_var(:,idx_q) = A(:,i).*vsVar;
                if idx_eq == i
                    coe_p_var(:,idx_p) = coe_p_var(:,idx_p) - DVS;
                    coe_q_var(:,idx_q) = coe_q_var(:,idx_q) + DVC;
                end
            end
        end

        function output = residue_system(System,tphi,TS_variable, ...
            TS_Dvariable,parameters,idx_eq)
            F = System{idx_eq};
            output = -TS_Dvariable(:,idx_eq) + tphi.*F(TS_variable,parameters);
        end

        function [coe_p_var,coe_q_var,coe_phi] = delta_coe_phi_var(p,q,phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
            coe_p_var = -(0:M).*[0,vs];
            coe_q_var = (1:M).*vc(2:end);
            coe_phi = -((0:M).^2.*vc*p+(1:M).^2.*vs*q);
        end

        function [coe_p_var,coe_q_var,coe_phi] = delta_coe_obs_phi(Derivative_obs, ...
                DDerivative_obs, p_var,q_var,phi,k)
            M = size(q_var,1);
            N = size(q_var,2);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
        
            coe_phi = 0;
            coe_p_var = zeros(size(p_var));
            coe_q_var = zeros(size(q_var));
        
            pt_var = vc*p_var+vs*q_var; % evaluated at a single time point
            pt_Dvar = -vs.*(1:M)*p_var(2:end,:) + vc(:,2:end).*(1:M)*q_var;
            for i = 1:N
                dxidphi = pt_Dvar(i);
                dykdxi = Derivative_obs(k,i).function(pt_var);
        
                coe_p_var(:,i) = coe_p_var(:,i) + (-dykdxi*(0:M).*[0 vs])';
                coe_q_var(:,i) = coe_q_var(:,i) + (dykdxi*(1:M).*vc(2:end))';
                coe_phi = coe_phi-dykdxi*((0:M).^2.*vc*p_var(:,i)+(1:M).^2.*vs*q_var(:,i));
                for j = 1:N
                    dxjdphi = pt_Dvar(j);
                    ddykdxixj = DDerivative_obs(k,i,j).function(pt_var);
                    coe_p_var(:,j) = coe_p_var(:,j) + ddykdxixj*dxidphi*vc';
                    coe_q_var(:,j) = coe_q_var(:,j) + ddykdxixj*dxidphi*vs';
                    coe_phi = coe_phi + dxidphi*ddykdxixj*dxjdphi;
                end
            end
        end

        function res = residue_phi_var(p,q,phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
            res = -vs.*(1:M)*p(2:end,:) + vc(:,2:end).*(1:M)*q;
        end

        function res = residue_phi_obs(Derivative_obs,p_var,q_var,phi,k)
            M = size(q_var,1);
            N = size(q_var,2);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
        
            res = 0;
            pt_var = vc*p_var+vs*q_var; % evaluated at a single time point
            pt_Dvar = -vs.*(1:M)*p_var(2:end,:) + vc(:,2:end).*(1:M)*q_var;
            for i = 1:N
                dykdxi = Derivative_obs(k,i).function(pt_var);
                dxidphi = pt_Dvar(i);
                res = res + dykdxi*dxidphi;
            end
        end

        function [coe_p_var,coe_q_var,coe_phi_var] = delta_coe_var_target(p,q,Phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(Phi,M,1);
        
            coe_p_var = vc;
            coe_q_var = vs;
        
            coe_phi_var = (1:M).*vc(2:end)*q-(0:M).*[0 vs]*p;
        end

        function [coe_p_var,coe_q_var] = delta_coe_observable(Derivative_obs,p_var,q_var,PV,L)
            assert(strcmpi(PV.name,'obs'))
            N = size(q_var,2);
            M = size(q_var,1);
            phi = (0:L-1)'*2*pi/L;

            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,L);

            TS_var = vc*p_var + vs*q_var;
        
            A = zeros(L,N);
            for i = 1:N
                temp = Derivative_obs(PV.idx,i).function(TS_var);
                A(:,i) = ones(L,1).*temp;
            end

            coe_p_var = zeros(L,N*(M+1));
            coe_q_var = zeros(L,N*M);
            for i = 1:N
                idx_p = (i-1)*(M+1) + (1:M+1);
                idx_q = (i-1)*M + (1:M);
                coe_p_var(:,idx_p) = A(:,i).*vc;
                coe_q_var(:,idx_q) = A(:,i).*vs;
            end
        end

        function output = var_target_curr(p,q,Phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(Phi,M,1);
            
            output = vc*p+vs*q;
        end

        function [coe_p_var,coe_q_var,coe_phi_obs] = ...
                 delta_coe_obs_target(Derivative_obs,p_var,q_var,phi,k)
            M = size(q_var,1);
            N = size(q_var,2);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
        
            coe_p_var = zeros(size(p_var));
            coe_q_var = zeros(size(q_var));
            coe_phi_obs = 0;
        
            pt_var = vc*p_var+vs*q_var;
            pt_Dvar = -vs.*(1:M)*p_var(2:end,:) + vc(:,2:end).*(1:M)*q_var;
            for i = 1:N
                dykdxi = Derivative_obs(k,i).function(pt_var);
                dxidphi = pt_Dvar(i);

                coe_p_var(:,i) = dykdxi*vc';
                coe_q_var(:,i) = dykdxi*vs';
                coe_phi_obs = coe_phi_obs + dykdxi*dxidphi;
            end
        end

        function output = obs_target_curr(obs,p_var,q_var,Phi,k)
            M = size(q_var,1);
            [vc,vs] = FMAM_ODE.Vec_CS(Phi,M,1);  
            pt_var = vc*p_var+vs*q_var;
        
            F = obs{k};
            output = F(pt_var);
        end

        function [coe_p_Psi,coe_q_Psi,coe_phi1,coe_phi2] = delta_coe_state_phase(p,q,Phi1,Phi2)
            M = size(q,1);
            [vc1,vs1] = FMAM_ODE.Vec_CS(Phi1,M,1);
            [vc2,vs2] = FMAM_ODE.Vec_CS(Phi2,M,1);
            coe_p_Psi = [Phi2-Phi1 (vs2-vs1).*((1:M).^(-1))];
            coe_q_Psi = -(vc2(:,2:end)-vc1(:,2:end)).*((1:M).^(-1));
            coe_phi1 = -(vc1*p+vs1*q);
            coe_phi2 = vc2*p+vs2*q;
        end
    end
end
