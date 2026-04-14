classdef fmam_state_ops
    % Canonical FMAM state helper boundary.
    % This is the only supported neutral numeric API for new production and test code.
    methods (Static)
        function TS_obs = getObs(obs,TS_var)
            TS_obs = zeros(size(TS_var,1),numel(obs));
            for k = 1:numel(obs)
                funcObs = obs{k};
                TS_obs(:,k) = funcObs(TS_var);
            end
        end

        function [phi_max,phi_min,amplitude,var_max,var_min] = FindExtreme(p,q,L,countMax,errMax)
            if nargin < 4 || isempty(countMax)
                countMax = fmam_state_defaults.countMax;
            end
            if nargin < 5 || isempty(errMax)
                errMax = fmam_state_defaults.errMax;
            end

            count = 0;
            while count < countMax
                phi = (0:L-1)'/L*2*pi;
                TS = fmam_state_ops.evaluateTrigSeries(phi,p,q);

                [var_max,index_max] = max(TS);
                [var_min,index_min] = min(TS);

                phi_max = phi(index_max);
                phi_min = phi(index_min);

                err_phiMax = abs(fmam_state_ops.residuePhiVar(p,q,phi_max));
                err_phiMin = abs(fmam_state_ops.residuePhiVar(p,q,phi_min));
                err = max(err_phiMax,err_phiMin);
                if err < errMax
                    break
                end

                L = L * 2;
                count = count + 1;
            end

            phi_max = phi(index_max);
            phi_min = phi(index_min);
            amplitude = (var_max - var_min) / 2;
        end

        function output = Trintegration(p,q,phi_L,phi_U)
            M = size(q,1);
            assert(isequal(size(phi_L),size(phi_U)))
            ptNum = numel(phi_L);
            invModes = (1:M).^(-1);
            [vc_U,vs_U] = fmam_state_ops.vecCS(phi_U,M,ptNum);
            [vc_L,vs_L] = fmam_state_ops.vecCS(phi_L,M,ptNum);
            U = p(1)*phi_U + (vs_U .* invModes) * p(2:end) - (vc_U(:,2:end) .* invModes) * q;
            L = p(1)*phi_L + (vs_L .* invModes) * p(2:end) - (vc_L(:,2:end) .* invModes) * q;
            output = U - L;
        end

        function [p,q] = projectFourierSeries(TS,M)
            L = size(TS,1);
            temp = fft(TS);
            p = 2*real(temp(1:M+1,:))/L;
            p(1,:) = p(1,:)/2;
            q = -2*imag(temp(2:M+1,:))/L;
        end

        function data = buildPrimaryBranchInterpolation(branchX, branchT, branchVar, branchObs)
            branchX = branchX(:);
            branchT = branchT(:);

            if branchX(1) > branchX(end)
                branchX = flipud(branchX);
                branchT = flipud(branchT);
                branchVar = flipud(branchVar);
                branchObs = flipud(branchObs);
            end

            [branchXUnique, uniqueIdx] = unique(branchX, 'stable');
            if numel(branchXUnique) < 2
                error('state:InvalidPrimaryWaveform', ...
                    'Primary branch interpolation requires at least two distinct waveform values.');
            end

            data = struct();
            data.x = branchXUnique;
            data.t = branchT(uniqueIdx);
            data.var = branchVar(uniqueIdx,:);
            data.obs = branchObs(uniqueIdx,:);
        end

        function [tNorm,TSVarNorm] = normalizePeriodicInputs(t,TS_var)
            % Canonical preprocessing for one-period trajectory samples.
            if size(TS_var,1) ~= numel(t)
                error('The number of time points must match the number of rows in TS_var.')
            end

            TSVarNorm = TS_var;
            if fmam_state_ops.hasRepeatedEndpoint(TSVarNorm)
                t(end) = [];
                TSVarNorm(end,:) = [];
            end
            tNorm = t - t(1);
        end

        function validateTrajectoryInputs(obs,params,t,TS_var,M,PV)
            % Canonical validation entry point for trajectory -> solverView reconstruction.
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end
            if ~iscell(obs)
                error('obs must be a cell array of observable functions.')
            end
            if ~isnumeric(TS_var) || ~ismatrix(TS_var) || isempty(TS_var)
                error('TS_var must be a non-empty numeric matrix.')
            end
            if any(~isfinite(TS_var(:)))
                error('TS_var must contain only finite values.')
            end
            if size(TS_var,1) ~= numel(t)
                error('The number of time points must match the number of rows in TS_var.')
            end
            if ~isnumeric(t) || isempty(t) || any(~isfinite(t))
                error('t must be a finite numeric vector.')
            end
            if any(diff(t) <= 0)
                error('t must be strictly increasing.')
            end
            if numel(t) < 3
                error('t must contain at least three samples over one period.')
            end
            if t(end) <= 0
                error('t must span one positive period after normalization.')
            end
            if ~isscalar(M) || M < 1 || M ~= floor(M)
                error('M must be a positive integer.')
            end
            if ~isnumeric(params) || any(~isfinite(params))
                error('Params must be a finite numeric vector.')
            end

            fmam_state_ops.validatePrimaryVariable(PV,size(TS_var,2),numel(obs));
        end

        function tf = hasRepeatedEndpoint(TS_var)
            if size(TS_var,1) < 2
                tf = false;
                return
            end

            scale = max(1,max(abs(TS_var(:))));
            tf = norm(TS_var(end,:) - TS_var(1,:),inf) <= 1e-8*scale;
        end

        function issue = timeMapInvariantIssue(PsiValues,p,q,numIntervals)
            issue = fmam_state_ops.positivePsiIssue(PsiValues);
            if ~issue.isValid
                return
            end

            issue = fmam_state_ops.positiveTimeIncrementIssue(p,q,numIntervals);
        end

        function issue = positivePsiIssue(PsiValues)
            issue = struct('isValid', true, 'identifier', '', 'message', '');
            if any(~isfinite(PsiValues))
                issue.isValid = false;
                issue.identifier = 'state:NonMonotoneTimeMap';
                issue.message = 'Psi(phi) must remain finite on the sampling grid.';
                return
            end

            tol = 1e-10*max(1,max(abs(PsiValues)));
            if any(PsiValues <= tol)
                issue.isValid = false;
                issue.identifier = 'state:NonMonotoneTimeMap';
                issue.message = 'Psi(phi) must stay strictly positive on the sampling grid.';
            end
        end

        function assertPositivePsi(PsiValues)
            issue = fmam_state_ops.positivePsiIssue(PsiValues);
            if ~issue.isValid
                error(issue.identifier,'%s',issue.message);
            end
        end

        function dt = timeIncrementsFromCoefficients(p,q,numIntervals)
            if nargin < 3 || isempty(numIntervals)
                numIntervals = fmam_state_defaults.LphiConst;
            end

            phiNodes = (0:numIntervals-1)'*2*pi/numIntervals;
            phiNext = [phiNodes(2:end); 2*pi];
            dt = fmam_state_ops.Trintegration(p,q,phiNodes,phiNext);
        end

        function issue = positiveTimeIncrementIssue(p,q,numIntervals)
            issue = struct('isValid', true, 'identifier', '', 'message', '');
            dt = fmam_state_ops.timeIncrementsFromCoefficients(p,q,numIntervals);
            if any(~isfinite(dt))
                issue.isValid = false;
                issue.identifier = 'state:NonMonotoneTimeMap';
                issue.message = 'The exact time increments must remain finite.';
                return
            end

            tol = 1e-10*max(1,max(abs(dt)));
            if any(dt <= tol)
                issue.isValid = false;
                issue.identifier = 'state:NonMonotoneTimeMap';
                issue.message = 'The exact time increments must stay strictly positive.';
            end
        end

        function assertPositiveTimeIncrements(p,q,numIntervals)
            issue = fmam_state_ops.positiveTimeIncrementIssue(p,q,numIntervals);
            if ~issue.isValid
                error(issue.identifier,'%s',issue.message);
            end
        end

        function values = evaluateTrigSeries(phi,p,q)
            [vc,vs] = fmam_state_ops.vecCS(phi,size(q,1),numel(phi));
            values = vc*p + vs*q;
        end

        function [p_variable,q_variable] = reprojectEqualTimeFourier(tSeries,TS_var,M)
            T = tSeries(end);
            L = size(TS_var,1);
            t_equal = (0:L-1)'/L*T;
            TS_var_equal = zeros(size(TS_var));

            for i = 1:size(TS_var,2)
                TS_var_equal(:,i) = spline(tSeries,TS_var(:,i),t_equal);
            end
            [p_variable,q_variable] = fmam_state_ops.projectFourierSeries(TS_var_equal,M);
        end

        function solverView = buildSolverViewFromTrajectory(obs,params,t,TS_var,M,PV,Lconst,Lphi,countMax,errMax)
            % Canonical trajectory -> solverView adapter used by constructors and fixtures.
            if nargin < 7 || isempty(Lconst)
                Lconst = fmam_state_defaults.Lconst;
            end
            if nargin < 8 || isempty(Lphi)
                Lphi = fmam_state_defaults.LphiConst;
            end
            if nargin < 9 || isempty(countMax)
                countMax = fmam_state_defaults.countMax;
            end
            if nargin < 10 || isempty(errMax)
                errMax = fmam_state_defaults.errMax;
            end

            t = t(:) - t(1);
            params = reshape(params,1,[]);
            TS_obs = fmam_state_ops.getObs(obs,TS_var);
            [~,~,p_Psi,q_Psi,p_var,q_var] = fmam_state_ops.reconstructSolverCoefficients( ...
                TS_var,TS_obs,t,M,PV,Lconst,Lphi);

            varDerived = fmam_state_ops.buildVariableDerivedState( ...
                p_var,q_var,p_Psi,q_Psi,Lphi,countMax,errMax);
            obsDerived = fmam_state_ops.buildObservableDerivedState( ...
                obs,p_var,q_var,p_Psi,q_Psi,Lphi,Lconst,countMax,errMax);

            solverView = struct( ...
                'params', params, ...
                'p_Psi', p_Psi, ...
                'q_Psi', q_Psi, ...
                'p_var', p_var, ...
                'q_var', q_var, ...
                'varPhiMax', varDerived.varPhiMax, ...
                'varPhiMin', varDerived.varPhiMin, ...
                'obsPhiMax', obsDerived.obsPhiMax, ...
                'obsPhiMin', obsDerived.obsPhiMin, ...
                'PV', PV);
            solverView = fmam_state_ops.normalizeSolverView(solverView);
        end

        function derived = buildDerivedView(obs,solverView,Lphi,Lconst,countMax,errMax)
            % Canonical derived view builder for downstream solver/postprocessing callers.
            if nargin < 3 || isempty(Lphi)
                Lphi = fmam_state_defaults.LphiConst;
            end
            if nargin < 4 || isempty(Lconst)
                Lconst = fmam_state_defaults.Lconst;
            end
            if nargin < 5 || isempty(countMax)
                countMax = fmam_state_defaults.countMax;
            end
            if nargin < 6 || isempty(errMax)
                errMax = fmam_state_defaults.errMax;
            end

            phi = (0:Lphi-1)'*2*pi/Lphi;
            derived = struct();
            derived.Psi = fmam_state_ops.evaluateTrigSeries(phi,solverView.p_Psi,solverView.q_Psi);
            dt = fmam_state_ops.timeIncrementsFromCoefficients( ...
                solverView.p_Psi,solverView.q_Psi,Lphi);
            derived.t = [0; cumsum(dt(1:end-1))];
            derived.TS_var = fmam_state_ops.evaluateTrigSeries(phi,solverView.p_var,solverView.q_var);
            derived.period = 2*pi*solverView.p_Psi(1);

            tScaled = derived.t / (derived.period / (2*pi));
            [derived.p_var_origin,derived.q_var_origin] = ...
                fmam_state_ops.reprojectEqualTimeFourier(tScaled,derived.TS_var,size(solverView.q_var,1));

            dimVar = size(solverView.p_var,2);
            derived.varMax = zeros(1,dimVar);
            derived.varMin = zeros(1,dimVar);
            for i = 1:dimVar
                derived.varMax(i) = fmam_state_ops.evaluateTrigSeries( ...
                    solverView.varPhiMax(i),solverView.p_var(:,i),solverView.q_var(:,i));
                derived.varMin(i) = fmam_state_ops.evaluateTrigSeries( ...
                    solverView.varPhiMin(i),solverView.p_var(:,i),solverView.q_var(:,i));
            end
            derived.varAmp = 0.5 * (derived.varMax - derived.varMin);
            if dimVar > 0
                tVarMax = fmam_state_ops.Trintegration( ...
                    solverView.p_Psi,solverView.q_Psi,zeros(dimVar,1),solverView.varPhiMax(:));
                derived.varPhase = repmat(tVarMax',dimVar,1) - repmat(tVarMax,1,dimVar);
            else
                derived.varPhase = zeros(0,0);
            end

            M = size(solverView.q_var,1);
            [derived.p_obs,derived.q_obs] = fmam_state_ops.buildObservableFourierCoefficients( ...
                obs,solverView.p_var,solverView.q_var,M,Lphi);
            if ~isempty(obs)
                derived.TS_obs = fmam_state_ops.evaluateTrigSeries(phi,derived.p_obs,derived.q_obs);
            else
                derived.TS_obs = zeros(size(derived.TS_var,1),0);
            end

            dimObs = numel(obs);
            derived.obsMax = zeros(1,dimObs);
            derived.obsMin = zeros(1,dimObs);
            for i = 1:dimObs
                derived.obsMax(i) = fmam_state_ops.observableTargetCurrent( ...
                    obs,solverView.p_var,solverView.q_var,solverView.obsPhiMax(i),i);
                derived.obsMin(i) = fmam_state_ops.observableTargetCurrent( ...
                    obs,solverView.p_var,solverView.q_var,solverView.obsPhiMin(i),i);
            end
            derived.obsAmp = 0.5 * (derived.obsMax - derived.obsMin);
            if dimObs > 0
                tObsMax = fmam_state_ops.Trintegration( ...
                    solverView.p_Psi,solverView.q_Psi,zeros(dimObs,1),solverView.obsPhiMax(:));
                derived.obsPhase = repmat(tObsMax',dimObs,1) - repmat(tObsMax,1,dimObs);
            else
                derived.obsPhase = zeros(0,0);
            end
        end

        function derived = buildVariableDerivedState(p_var,q_var,p_Psi,q_Psi,Lphi,countMax,errMax)
            dimSys = size(p_var,2);
            Phi_max = zeros(1,dimSys);
            Phi_min = zeros(1,dimSys);
            Amp = zeros(1,dimSys);
            VariableMax = zeros(1,dimSys);
            VariableMin = zeros(1,dimSys);

            for i = 1:dimSys
                [Phi_max(i),Phi_min(i),Amp(i),VariableMax(i),VariableMin(i)] = ...
                    fmam_state_ops.FindExtreme(p_var(:,i),q_var(:,i),Lphi,countMax,errMax);
            end
            tMax = fmam_state_ops.Trintegration(p_Psi,q_Psi,zeros(dimSys,1),Phi_max');
            Phase = repmat(tMax',dimSys,1) - repmat(tMax,1,dimSys);

            derived = struct();
            derived.varPhiMax = Phi_max;
            derived.varPhiMin = Phi_min;
            derived.varAmp = Amp;
            derived.varMax = VariableMax;
            derived.varMin = VariableMin;
            derived.varPhase = Phase;
        end

        function [p_obs,q_obs] = buildObservableFourierCoefficients(obs,p_var,q_var,M,Lphi)
            if isempty(obs)
                p_obs = zeros(M+1,0);
                q_obs = zeros(M,0);
                return
            end

            phi = (0:Lphi-1)'*2*pi/Lphi;
            TS_var = fmam_state_ops.evaluateTrigSeries(phi,p_var,q_var);
            observableSeries = fmam_state_ops.getObs(obs,TS_var);
            [p_obs,q_obs] = fmam_state_ops.projectFourierSeries(observableSeries,M);
        end

        function derived = buildObservableDerivedState(obs,p_var,q_var,p_Psi,q_Psi,Lphi,Lconst,countMax,errMax)
            M = size(q_var,1);
            [p_obs,q_obs] = fmam_state_ops.buildObservableFourierCoefficients(obs,p_var,q_var,M,Lphi);
            dimObs = size(p_obs,2);

            Phi_max = zeros(1,dimObs);
            Phi_min = zeros(1,dimObs);
            Amp = zeros(1,dimObs);
            ObservableMax = zeros(1,dimObs);
            ObservableMin = zeros(1,dimObs);

            for i = 1:dimObs
                [Phi_max(i),Phi_min(i),Amp(i),ObservableMax(i),ObservableMin(i)] = ...
                    fmam_state_ops.FindExtreme(p_obs(:,i),q_obs(:,i),Lconst,countMax,errMax);
            end
            tMax = fmam_state_ops.Trintegration(p_Psi,q_Psi,zeros(dimObs,1),Phi_max');
            Phase = repmat(tMax',dimObs,1) - repmat(tMax,1,dimObs);

            derived = struct();
            derived.p_obs = p_obs;
            derived.q_obs = q_obs;
            derived.obsPhiMax = Phi_max;
            derived.obsPhiMin = Phi_min;
            derived.obsAmp = Amp;
            derived.obsMax = ObservableMax;
            derived.obsMin = ObservableMin;
            derived.obsPhase = Phase;
        end

        function [a,b,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = ...
                reconstructSolverCoefficients(TS_variable,TS_observable,t,M,PV,Lconst,Lphi)
            dim = size(TS_variable,2);
            n = size(TS_observable,2);
            fmam_state_ops.validatePrimaryVariable(PV,dim,n);
            p_observable = [];
            q_observable = [];

            T = t(end);
            if T <= 0
                error('state:InvalidTimeRange', ...
                    't must span one positive period.')
            end

            phi = (0:Lphi-1)'*2*pi/Lphi;
            t1 = linspace(0,T,Lconst+1)';
            t1(end) = [];
            TS_variable_1 = zeros(Lconst,dim);
            TS_observable_1 = zeros(Lconst,n);
            for i = 1:dim
                TS_variable_1(:,i) = spline(t,TS_variable(:,i),t1);
            end

            for i = 1:n
                TS_observable_1(:,i) = spline(t,TS_observable(:,i),t1);
            end

            if strcmpi(PV.name,'var')
                X = TS_variable_1(:,PV.idx);
            elseif strcmpi(PV.name,'obs')
                X = TS_observable_1(:,PV.idx);
            else
                error('Please check the class of the parimary variable.')
            end

            [~,I1] = max(X);
            TS_variable_1 = [TS_variable_1(I1:end,:);TS_variable_1(1:I1-1,:)];
            X = [X(I1:end);X(1:I1-1)];
            if n > 0
                TS_observable_1 = [TS_observable_1(I1:end,:);TS_observable_1(1:I1-1,:)];
            end
            t1 = [t1(I1:end)-t1(I1);t1(1:I1-1)+T-t1(I1)];

            [~,I2] = min(X);
            a = (max(X)-min(X))/2;
            b = (max(X)+min(X))/2;
            if a <= 0
                error('state:InvalidPrimaryWaveform', ...
                    'Primary variable must have a nonzero amplitude over one period.');
            end

            leftBranch = 1:I2;
            rightBranch = I2:numel(X);
            if numel(leftBranch) < 2 || numel(rightBranch) < 2
                error('state:InvalidPrimaryWaveform', ...
                    ['Primary variable must admit usable max-to-min and min-to-max ', ...
                     'branches over one period.']);
            end

            Atrans = zeros(Lphi,Lphi+1);
            Atrans(1,1) = -3;Atrans(1,2) = 4;Atrans(1,3) = -1;
            for j = 2:Lphi
                Atrans(j,j+1) = 1;
                Atrans(j,j-1) = -1;
            end
            Atrans = Atrans/(4*pi/Lphi);

            t_phi = zeros(Lphi,1);
            TS_variable_phi = zeros(Lphi,dim);
            TS_observable_phi = zeros(Lphi,n);
            leftData = fmam_state_ops.buildPrimaryBranchInterpolation( ...
                X(leftBranch),t1(leftBranch),TS_variable_1(leftBranch,:),TS_observable_1(leftBranch,:));
            rightData = fmam_state_ops.buildPrimaryBranchInterpolation( ...
                X(rightBranch),t1(rightBranch),TS_variable_1(rightBranch,:),TS_observable_1(rightBranch,:));
            for i = 1:Lphi
                targetValue = a*cos(phi(i))+b;
                if i <= round(Lphi/2)
                    branchData = leftData;
                else
                    branchData = rightData;
                end

                t_phi(i) = interp1(branchData.x,branchData.t,targetValue,'pchip');
                TS_variable_phi(i,:) = interp1(branchData.x,branchData.var,targetValue,'pchip');
                if n > 0
                    TS_observable_phi(i,:) = interp1(branchData.x,branchData.obs,targetValue,'pchip');
                end
            end

            t_phi = [t_phi;T];
            Psi = Atrans*t_phi;

            [p_Psi,q_Psi] = fmam_state_ops.projectFourierSeries(Psi,M - 1);
            fmam_state_ops.assertPositivePsi(fmam_state_ops.evaluateTrigSeries(phi,p_Psi,q_Psi));
            fmam_state_ops.assertPositiveTimeIncrements(p_Psi,q_Psi,Lphi);

            [p_variable,q_variable] = fmam_state_ops.projectFourierSeries(TS_variable_phi,M);
            if n > 0
                [p_observable,q_observable] = fmam_state_ops.projectFourierSeries(TS_observable_phi,M);
            end
            if strcmpi(PV.name,'var')
                p_variable(1,PV.idx) = b;
                p_variable(2,PV.idx) = a;
                p_variable(3:end,PV.idx) = zeros(M-1,1);
                q_variable(:,PV.idx) = zeros(M,1);
            elseif strcmpi(PV.name,'obs')
                p_observable(3:end,PV.idx) = zeros(M-1,1);
                q_observable(:,PV.idx) = zeros(M,1);
            end
        end

        function validatePrimaryVariable(PV,dimVar,dimObs)
            if ~isstruct(PV) || ~all(isfield(PV,{'name','idx'}))
                error('PV must be a struct with fields ''name'' and ''idx''.')
            end
            if ~ischar(PV.name) && ~(isstring(PV.name) && isscalar(PV.name))
                error('PV.name must be ''var'' or ''obs''.')
            end
            if ~isscalar(PV.idx) || PV.idx < 1 || PV.idx ~= floor(PV.idx)
                error('PV.idx must be a positive integer.')
            end

            switch lower(char(PV.name))
                case 'var'
                    if PV.idx > dimVar
                        error('PV.idx exceeds the number of state variables.')
                    end
                case 'obs'
                    if dimObs == 0
                        error('PV.name cannot be ''obs'' when no observables are provided.')
                    end
                    if PV.idx > dimObs
                        error('PV.idx exceeds the number of observables.')
                    end
                otherwise
                    error('Please check the class of the parimary variable.')
            end
        end

        function normalized = normalizeSolverView(view,dimVar,dimObs,dimParams,truncationOrder)
            requiredFields = {'params','p_Psi','q_Psi','p_var','q_var', ...
                'varPhiMax','varPhiMin','obsPhiMax','obsPhiMin','PV'};
            if ~isstruct(view) || ~all(isfield(view,requiredFields))
                error('FMAM_ODE:InvalidSolverView', ...
                    'solverView must contain the fields %s.', strjoin(requiredFields, ', '))
            end

            if nargin < 5 || isempty(truncationOrder)
                if isfield(view,'truncationOrder') && ~isempty(view.truncationOrder)
                    truncationOrder = double(view.truncationOrder);
                elseif isfield(view,'p_var') && ~isempty(view.p_var)
                    truncationOrder = size(view.p_var,1) - 1;
                elseif isfield(view,'q_var') && ~isempty(view.q_var)
                    truncationOrder = size(view.q_var,1);
                else
                    truncationOrder = size(view.p_var,1) - 1;
                end
            end
            if nargin < 2 || isempty(dimVar)
                dimVar = size(view.p_var,2);
            end
            if nargin < 3 || isempty(dimObs)
                dimObs = size(reshape(view.obsPhiMax,1,[]),2);
            end
            if nargin < 4 || isempty(dimParams)
                dimParams = numel(view.params);
            end

            normalized = struct();
            normalized.params = reshape(view.params,1,[]);
            normalized.p_Psi = view.p_Psi;
            normalized.q_Psi = view.q_Psi;
            normalized.p_var = view.p_var;
            normalized.q_var = view.q_var;
            normalized.varPhiMax = reshape(view.varPhiMax,1,[]);
            normalized.varPhiMin = reshape(view.varPhiMin,1,[]);
            normalized.obsPhiMax = reshape(view.obsPhiMax,1,[]);
            normalized.obsPhiMin = reshape(view.obsPhiMin,1,[]);
            normalized.PV = view.PV;
            fmam_state_ops.validatePrimaryVariable(normalized.PV,dimVar,dimObs);
            normalized.propertySizes = fmam_state_ops.solverPropertySizesFromArrays( ...
                normalized.params,normalized.p_Psi,normalized.q_Psi, ...
                normalized.p_var,normalized.q_var,normalized.varPhiMax, ...
                normalized.varPhiMin,normalized.obsPhiMax,normalized.obsPhiMin);
            normalized.dimSys = dimVar;
            normalized.dimObs = dimObs;
            normalized.dimParams = dimParams;
            normalized.truncationOrder = truncationOrder;

            if normalized.dimParams ~= numel(normalized.params)
                error('FMAM_ODE:InvalidSolverView', ...
                    'solverView.params must match the configured parameter dimension.')
            end
            if ~isequal(size(normalized.p_Psi),[normalized.truncationOrder, 1]) || ...
                    ~isequal(size(normalized.q_Psi),[max(normalized.truncationOrder - 1, 0), 1])
                error('FMAM_ODE:InvalidSolverView', ...
                    'solverView p_Psi/q_Psi shapes must match truncationOrder.')
            end
            if ~isequal(size(normalized.p_var),[normalized.truncationOrder + 1, normalized.dimSys]) || ...
                    ~isequal(size(normalized.q_var),[normalized.truncationOrder, normalized.dimSys])
                error('FMAM_ODE:InvalidSolverView', ...
                    'solverView p_var/q_var shapes must match the configured system dimension.')
            end
            if ~isequal(size(normalized.varPhiMax),[1, normalized.dimSys]) || ...
                    ~isequal(size(normalized.varPhiMin),[1, normalized.dimSys])
                error('FMAM_ODE:InvalidSolverView', ...
                    'solverView variable phase arrays must be row vectors of length dimSys.')
            end
            if ~isequal(size(normalized.obsPhiMax),[1, normalized.dimObs]) || ...
                    ~isequal(size(normalized.obsPhiMin),[1, normalized.dimObs])
                error('FMAM_ODE:InvalidSolverView', ...
                    'solverView observable phase arrays must be row vectors of length dimObs.')
            end
        end

        function sizes = solverPropertySizesFromArrays(params,p_Psi,q_Psi,p_var,q_var,varPhiMax,varPhiMin,obsPhiMax,obsPhiMin)
            sizes = struct( ...
                'params', size(params), ...
                'p_Psi', size(p_Psi), ...
                'q_Psi', size(q_Psi), ...
                'p_var', size(p_var), ...
                'q_var', size(q_var), ...
                'varPhiMax', size(varPhiMax), ...
                'varPhiMin', size(varPhiMin), ...
                'obsPhiMax', size(obsPhiMax), ...
                'obsPhiMin', size(obsPhiMin));
        end

        function solverView = solverViewFromState(stateObj)
            % Explicit compat adapter from rich state objects into canonical solverView structs.
            solverView = fmam_state_ops.normalizeSolverView(struct( ...
                'params', reshape(stateObj.params,1,[]), ...
                'p_Psi', stateObj.p_Psi, ...
                'q_Psi', stateObj.q_Psi, ...
                'p_var', stateObj.p_var, ...
                'q_var', stateObj.q_var, ...
                'varPhiMax', reshape(stateObj.varPhiMax,1,[]), ...
                'varPhiMin', reshape(stateObj.varPhiMin,1,[]), ...
                'obsPhiMax', reshape(stateObj.obsPhiMax,1,[]), ...
                'obsPhiMin', reshape(stateObj.obsPhiMin,1,[]), ...
                'PV', stateObj.PV, ...
                'truncationOrder', stateObj.truncationOrder), ...
                stateObj.dimSys,stateObj.dimObs,stateObj.dimParams,stateObj.truncationOrder);
        end

        function derived = derivedViewFromState(stateObj)
            % Explicit compat adapter from rich state objects into canonical derived structs.
            derived = struct( ...
                'Psi', stateObj.Psi, ...
                't', stateObj.t, ...
                'TS_var', stateObj.TS_var, ...
                'TS_obs', stateObj.TS_obs, ...
                'period', stateObj.period, ...
                'p_var_origin', stateObj.p_var_origin, ...
                'q_var_origin', stateObj.q_var_origin, ...
                'p_obs', stateObj.p_obs, ...
                'q_obs', stateObj.q_obs, ...
                'varAmp', stateObj.varAmp, ...
                'obsAmp', stateObj.obsAmp, ...
                'varMax', stateObj.varMax, ...
                'obsMax', stateObj.obsMax, ...
                'varMin', stateObj.varMin, ...
                'obsMin', stateObj.obsMin, ...
                'varPhase', stateObj.varPhase, ...
                'obsPhase', stateObj.obsPhase);
        end

        function solverView = solverViewFromSnapshot(snapshot)
            % Canonical extraction of solver coefficients from legacy rich snapshots.
            solverView = fmam_state_ops.normalizeSolverView(struct( ...
                'params', snapshot.params, ...
                'p_Psi', snapshot.p_Psi, ...
                'q_Psi', snapshot.q_Psi, ...
                'p_var', snapshot.p_var, ...
                'q_var', snapshot.q_var, ...
                'varPhiMax', snapshot.varPhiMax, ...
                'varPhiMin', snapshot.varPhiMin, ...
                'obsPhiMax', snapshot.obsPhiMax, ...
                'obsPhiMin', snapshot.obsPhiMin, ...
                'PV', snapshot.PV));
        end

        function derived = derivedViewFromSnapshot(snapshot,Lphi)
            % Canonical extraction of derived arrays from legacy rich snapshots.
            if nargin < 2 || isempty(Lphi)
                Lphi = fmam_state_defaults.LphiConst;
            end
            phi = (0:Lphi-1)'*2*pi/Lphi;
            Psi = fmam_state_ops.evaluateTrigSeries(phi,snapshot.p_Psi,snapshot.q_Psi);
            dt = fmam_state_ops.timeIncrementsFromCoefficients(snapshot.p_Psi,snapshot.q_Psi,Lphi);
            TS_var = fmam_state_ops.evaluateTrigSeries(phi,snapshot.p_var,snapshot.q_var);
            derived = struct();
            derived.Psi = Psi;
            derived.t = [0; cumsum(dt(1:end-1))];
            derived.TS_var = TS_var;
            derived.TS_obs = fmam_state_ops.evaluateTrigSeries(phi,snapshot.p_obs,snapshot.q_obs);
            derived.period = snapshot.period;
            derived.p_var_origin = snapshot.p_var_origin;
            derived.q_var_origin = snapshot.q_var_origin;
            derived.p_obs = snapshot.p_obs;
            derived.q_obs = snapshot.q_obs;
            derived.varAmp = snapshot.varAmp;
            derived.obsAmp = snapshot.obsAmp;
            derived.varMax = snapshot.varMax;
            derived.obsMax = snapshot.obsMax;
            derived.varMin = snapshot.varMin;
            derived.obsMin = snapshot.obsMin;
            derived.varPhase = snapshot.varPhase;
            derived.obsPhase = snapshot.obsPhase;
        end

        function snapshot = buildStateSnapshotFromViews(solverView,derived)
            % Canonical view -> rich snapshot adapter for legacy state restore/rebuild paths.
            [a,b] = fmam_state_ops.primaryAmplitudeAndCenter(derived,solverView.PV);
            snapshot = struct( ...
                'PV', solverView.PV, ...
                'params', reshape(solverView.params,1,[]), ...
                'p_Psi', solverView.p_Psi, ...
                'q_Psi', solverView.q_Psi, ...
                'p_var', solverView.p_var, ...
                'q_var', solverView.q_var, ...
                'p_var_origin', derived.p_var_origin, ...
                'q_var_origin', derived.q_var_origin, ...
                'p_obs', derived.p_obs, ...
                'q_obs', derived.q_obs, ...
                'a', a, ...
                'b', b, ...
                'varPhiMax', reshape(solverView.varPhiMax,1,[]), ...
                'obsPhiMax', reshape(solverView.obsPhiMax,1,[]), ...
                'varPhiMin', reshape(solverView.varPhiMin,1,[]), ...
                'obsPhiMin', reshape(solverView.obsPhiMin,1,[]), ...
                'varAmp', derived.varAmp, ...
                'obsAmp', derived.obsAmp, ...
                'period', derived.period, ...
                'varMax', derived.varMax, ...
                'obsMax', derived.obsMax, ...
                'varMin', derived.varMin, ...
                'obsMin', derived.obsMin, ...
                'varPhase', derived.varPhase, ...
                'obsPhase', derived.obsPhase);
        end

        function [a,b] = primaryAmplitudeAndCenter(derived,PV)
            if strcmpi(PV.name,'var')
                a = derived.varAmp(PV.idx);
                b = 0.5 * (derived.varMax(PV.idx) + derived.varMin(PV.idx));
            else
                a = derived.obsAmp(PV.idx);
                b = 0.5 * (derived.obsMax(PV.idx) + derived.obsMin(PV.idx));
            end
        end
    end

    methods (Static, Access = private)
        function value = observableTargetCurrent(obs,p_var,q_var,Phi,k)
            pt_var = fmam_state_ops.evaluateTrigSeries(Phi,p_var,q_var);
            value = obs{k}(pt_var);
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
            [vc,vs] = fmam_state_ops.vecCS(phi,M,1);
            res = -vs .* (1:M) * p(2:end,:) + vc(:,2:end) .* (1:M) * q;
        end
    end
end
