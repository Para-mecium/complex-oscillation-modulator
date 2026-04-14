classdef state < handle
    properties (Constant)
        Lconst = 20000
        LphiConst = 500
        countMax = 20
        errMax = 1e-2
    end

    properties
        PV
        checkPsiNonnegative = true
    end

    properties (Access = private)
        timeMapWarningSignature = ''
    end

    properties (SetAccess=private)
        obs
        % category 1
        params
        p_Psi
        q_Psi
        p_var
        q_var
        p_var_origin
        q_var_origin
        p_obs
        q_obs
        % category 2
        a
        b
        varPhiMax = []
        obsPhiMax = []
        varPhiMin = []
        obsPhiMin = []
        varAmp = []
        obsAmp = []
        period = []
        varMax = []
        obsMax = []
        varMin = []
        obsMin = []
        varPhase = []
        obsPhase = []
    end

    properties(Dependent)
        dimSys
        dimObs
        dimParams
        truncationOrder
        Atrans

        t
        phi
        Psi
        TS_var
        TS_obs
    end
    methods
        function obj = state(obs,Params,t,TS_var,M,PV)
            if nargin == 0
                return
            end
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end
            t = t(:);
            Params = Params(:)';
            [t,TS_var] = state.normalizePeriodicInputs(t,TS_var);
            state.validateInputs(obs,Params,t,TS_var,M,PV);

            obj.obs = obs;
            solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
                obs,Params,t,TS_var,M,PV,obj.Lconst,obj.LphiConst,obj.countMax,obj.errMax);
            obj.params = solverView.params;
            obj.p_Psi = solverView.p_Psi;
            obj.q_Psi = solverView.q_Psi;
            obj.p_var = solverView.p_var;
            obj.q_var = solverView.q_var;
            obj.varPhiMax = solverView.varPhiMax;
            obj.varPhiMin = solverView.varPhiMin;
            obj.obsPhiMax = solverView.obsPhiMax;
            obj.obsPhiMin = solverView.obsPhiMin;
            obj.PV = solverView.PV;
            obj.refreshDerivedState();
        end

        function assertTimeMapInvariant(obj)
            issue = fmam_state_ops.timeMapInvariantIssue(obj.Psi,obj.p_Psi,obj.q_Psi,obj.LphiConst);
            obj.warnOnTimeMapIssue(issue);
        end
        
        function prop = get.dimSys(obj)
            prop = size(obj.p_var,2);
        end

        function prop = get.dimObs(obj)
            prop = size(obj.p_obs,2);
        end

        function prop = get.dimParams(obj)
            prop = size(obj.params,2);
        end

        function prop = get.truncationOrder(obj)
            prop = size(obj.q_var,1);
        end

        function prop = get.TS_var(obj)
            prop = fmam_state_ops.evaluateTrigSeries(obj.phi,obj.p_var,obj.q_var);
        end

        function prop = get.TS_obs(obj)
            prop = fmam_state_ops.getObs(obj.obs,obj.TS_var);
        end

        function prop = get.Psi(obj)
            prop = fmam_state_ops.evaluateTrigSeries(obj.phi,obj.p_Psi,obj.q_Psi);
        end

        function prop = get.phi(obj)
            L = obj.LphiConst;
            prop = (0:L-1)'*2*pi/L;
        end

        function prop = get.Atrans(obj)
            L = obj.LphiConst;
            prop = zeros(L,L+1);
            prop(1,1) = -3;prop(1,2) = 4;prop(1,3) = -1;
            for j = 2:L
                prop(j,j+1) = 1;
                prop(j,j-1) = -1;
            end
            prop = prop/(4*pi/L);
        end

        function prop = get.t(obj)
            issue = fmam_state_ops.timeMapInvariantIssue(obj.Psi,obj.p_Psi,obj.q_Psi,obj.LphiConst);
            obj.warnOnTimeMapIssue(issue);
            dt = fmam_state_ops.timeIncrementsFromCoefficients(obj.p_Psi,obj.q_Psi,obj.LphiConst);
            prop = [0; cumsum(dt(1:end-1))];
        end

        function updateVar2(obj)
            derived = fmam_state_ops.buildVariableDerivedState( ...
                obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi, ...
                obj.LphiConst,obj.countMax,obj.errMax);

            obj.varPhiMax = derived.varPhiMax;
            obj.varPhiMin = derived.varPhiMin;
            obj.varAmp = derived.varAmp;
            obj.varMax = derived.varMax;
            obj.varMin = derived.varMin;
            obj.varPhase = derived.varPhase;
        end

        function updateObs2(obj)
            if obj.dimObs == 0
                return
            end
            derived = fmam_state_ops.buildObservableDerivedState( ...
                obj.obs,obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi, ...
                obj.LphiConst,obj.Lconst,obj.countMax,obj.errMax);

            obj.p_obs = derived.p_obs;
            obj.q_obs = derived.q_obs;
            obj.obsPhiMax = derived.obsPhiMax;
            obj.obsPhiMin = derived.obsPhiMin;
            obj.obsAmp = derived.obsAmp;
            obj.obsMax = derived.obsMax;
            obj.obsMin = derived.obsMin;
            obj.obsPhase = derived.obsPhase;
        end

        function updatePeriod(obj)
            obj.period = 2*pi*obj.p_Psi(1);
        end

        function add(obj,propname,idx,val)
            if strcmpi(idx,'all')
                if ~isequal(size(obj.(propname)),size(val))
                    val = val';
                end
                obj.(propname) = obj.(propname) + val;
            else
                obj.(propname)(idx) = obj.(propname)(idx) + val;
            end
        end

        function minus(obj,propname,idx,val)
            if strcmpi(idx,'all')
                if ~isequal(size(obj.(propname)),size(val))
                    val = val';
                end
                obj.(propname) = obj.(propname) - val;
            else
                obj.(propname)(idx) = obj.(propname)(idx) - val;
            end
        end

        function changeTruncOrder(obj,M)
            M1 = obj.truncationOrder;
            if M1 >= M
                obj.p_var = obj.p_var(1:M+1,:);
                obj.q_var = obj.q_var(1:M,:);
                obj.p_Psi = obj.p_Psi(1:M,:);
                obj.q_Psi = obj.q_Psi(1:max(M-1,0),:);
            else
                obj.p_var = [obj.p_var;zeros(M-M1,obj.dimSys)];
                obj.q_var = [obj.q_var;zeros(M-M1,obj.dimSys)];
                obj.p_Psi = [obj.p_Psi;zeros(M-M1,1)];
                obj.q_Psi = [obj.q_Psi;zeros(M-M1,1)];
            end
            obj.refreshDerivedState();
        end

        function TSplot(obj,propname,idx,timedomain)
            plot(obj.(timedomain),[obj.(propname)(:,idx)])

            xlabel(timedomain)
        end

        function updatePV(obj,PV)
            state.validatePrimaryVariable(PV,obj.dimSys,size(obj.obs,2));
            solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
                obj.obs,obj.params,obj.t,obj.TS_var,obj.truncationOrder,PV, ...
                obj.Lconst,obj.LphiConst,obj.countMax,obj.errMax);
            obj.PV = solverView.PV;
            obj.params = solverView.params;
            obj.p_Psi = solverView.p_Psi;
            obj.q_Psi = solverView.q_Psi;
            obj.p_var = solverView.p_var;
            obj.q_var = solverView.q_var;
            obj.varPhiMax = solverView.varPhiMax;
            obj.varPhiMin = solverView.varPhiMin;
            obj.obsPhiMax = solverView.obsPhiMax;
            obj.obsPhiMin = solverView.obsPhiMin;
            obj.refreshDerivedState();
        end

        function refreshDerivedState(obj)
            obj.assertTimeMapInvariant();
            varDerived = fmam_state_ops.buildVariableDerivedState( ...
                obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi, ...
                obj.LphiConst,obj.countMax,obj.errMax);
            obsDerived = fmam_state_ops.buildObservableDerivedState( ...
                obj.obs,obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi, ...
                obj.LphiConst,obj.Lconst,obj.countMax,obj.errMax);

            solverView = struct( ...
                'params', obj.params, ...
                'p_Psi', obj.p_Psi, ...
                'q_Psi', obj.q_Psi, ...
                'p_var', obj.p_var, ...
                'q_var', obj.q_var, ...
                'varPhiMax', varDerived.varPhiMax, ...
                'varPhiMin', varDerived.varPhiMin, ...
                'obsPhiMax', obsDerived.obsPhiMax, ...
                'obsPhiMin', obsDerived.obsPhiMin, ...
                'PV', obj.PV);
            derived = fmam_state_ops.buildDerivedView( ...
                obj.obs,solverView,obj.LphiConst,obj.Lconst,obj.countMax,obj.errMax);

            obj.varPhiMax = solverView.varPhiMax;
            obj.varPhiMin = solverView.varPhiMin;
            obj.obsPhiMax = solverView.obsPhiMax;
            obj.obsPhiMin = solverView.obsPhiMin;
            obj.p_var_origin = derived.p_var_origin;
            obj.q_var_origin = derived.q_var_origin;
            obj.p_obs = derived.p_obs;
            obj.q_obs = derived.q_obs;
            obj.varAmp = derived.varAmp;
            obj.obsAmp = derived.obsAmp;
            obj.period = derived.period;
            obj.varMax = derived.varMax;
            obj.obsMax = derived.obsMax;
            obj.varMin = derived.varMin;
            obj.obsMin = derived.obsMin;
            obj.varPhase = derived.varPhase;
            obj.obsPhase = derived.obsPhase;
            [obj.a,obj.b] = state.primaryAmplitudeAndCenter(derived,obj.PV);
        end

        function updateFCOrigin(obj)
            derived = fmam_state_ops.buildDerivedView( ...
                obj.obs,obj.currentSolverView(),obj.LphiConst,obj.Lconst,obj.countMax,obj.errMax);
            obj.p_var_origin = derived.p_var_origin;
            obj.q_var_origin = derived.q_var_origin;
        end
    end

    methods (Access = private)
        function loadSolverSnapshot(obj,snapshot)
            fieldNames = state.solverStateFields();
            if ~isstruct(snapshot) || ~all(isfield(snapshot,fieldNames))
                error('state:InvalidSnapshot', ...
                    'snapshot must contain the full solver state.');
            end

            for i = 1:numel(fieldNames)
                fieldName = fieldNames{i};
                obj.(fieldName) = snapshot.(fieldName);
            end
        end

        function warnOnTimeMapIssue(obj,issue)
            if ~obj.checkPsiNonnegative
                obj.timeMapWarningSignature = '';
                return
            end
            if issue.isValid
                obj.timeMapWarningSignature = '';
                return
            end

            signature = [issue.identifier '|' issue.message];
            if ~strcmp(obj.timeMapWarningSignature,signature)
                warning(issue.identifier,'%s',issue.message);
                obj.timeMapWarningSignature = signature;
            end
        end

        function refreshObservableFourierCoefficients(obj)
            [pObs,qObs] = fmam_state_ops.buildObservableFourierCoefficients( ...
                obj.obs,obj.p_var,obj.q_var,obj.truncationOrder,obj.LphiConst);
            obj.p_obs = pObs;
            obj.q_obs = qObs;
        end

        function solverView = currentSolverView(obj)
            solverView = struct( ...
                'params', obj.params, ...
                'p_Psi', obj.p_Psi, ...
                'q_Psi', obj.q_Psi, ...
                'p_var', obj.p_var, ...
                'q_var', obj.q_var, ...
                'varPhiMax', obj.varPhiMax, ...
                'varPhiMin', obj.varPhiMin, ...
                'obsPhiMax', obj.obsPhiMax, ...
                'obsPhiMin', obj.obsPhiMin, ...
                'PV', obj.PV);
        end
    end

    
    methods(Static)
        function TS_obs = getObs(obs,TS_var)
            TS_obs = fmam_state_ops.getObs(obs,TS_var);
        end
        function [phi_max,phi_min,amplitude,var_max,var_min] = FindExtreme(p,q,L)
            [phi_max,phi_min,amplitude,var_max,var_min] = ...
                fmam_state_ops.FindExtreme(p,q,L,state.countMax,state.errMax);
        end

        function output = Trintegration(p,q,phi_L,phi_U)
            output = fmam_state_ops.Trintegration(p,q,phi_L,phi_U);
        end

        function [a,b,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = ...
            fourierCoeffs(M,t,TS_variable,TS_observable,PV)
            [a,b,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = ...
                fmam_state_ops.reconstructSolverCoefficients( ...
                    TS_variable,TS_observable,t(:)-t(1),M,PV,state.Lconst,state.LphiConst);
        end

        function [p,q] = projectFourierSeries(TS,M)
            [p,q] = fmam_state_ops.projectFourierSeries(TS,M);
        end

        function data = buildPrimaryBranchInterpolation(branchX, branchT, branchVar, branchObs)
            data = fmam_state_ops.buildPrimaryBranchInterpolation(branchX, branchT, branchVar, branchObs);
        end

        function validateInputs(obs,Params,t,TS_var,M,PV)
            fmam_state_ops.validateTrajectoryInputs(obs,Params,t,TS_var,M,PV);
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

        function [tNorm,TSVarNorm] = normalizePeriodicInputs(t,TS_var)
            [tNorm,TSVarNorm] = fmam_state_ops.normalizePeriodicInputs(t,TS_var);
        end

        function tf = hasRepeatedEndpoint(TS_var)
            tf = fmam_state_ops.hasRepeatedEndpoint(TS_var);
        end

        function issue = timeMapInvariantIssue(PsiValues,p,q,numIntervals)
            issue = fmam_state_ops.timeMapInvariantIssue(PsiValues,p,q,numIntervals);
        end

        function issue = positivePsiIssue(PsiValues)
            issue = fmam_state_ops.positivePsiIssue(PsiValues);
        end

        function assertPositivePsi(PsiValues)
            fmam_state_ops.assertPositivePsi(PsiValues);
        end

        function dt = timeIncrementsFromCoefficients(p,q,numIntervals)
            if nargin < 3 || isempty(numIntervals)
                numIntervals = state.LphiConst;
            end
            dt = fmam_state_ops.timeIncrementsFromCoefficients(p,q,numIntervals);
        end

        function issue = positiveTimeIncrementIssue(p,q,numIntervals)
            issue = fmam_state_ops.positiveTimeIncrementIssue(p,q,numIntervals);
        end

        function assertPositiveTimeIncrements(p,q,numIntervals)
            fmam_state_ops.assertPositiveTimeIncrements(p,q,numIntervals);
        end

        function values = evaluateTrigSeries(phi,p,q)
            values = fmam_state_ops.evaluateTrigSeries(phi,p,q);
        end

        function fieldNames = solverStateFields()
            fieldNames = { ...
                'PV','params','p_Psi','q_Psi','p_var','q_var', ...
                'p_var_origin','q_var_origin','p_obs','q_obs', ...
                'a','b','varPhiMax','obsPhiMax','varPhiMin','obsPhiMin', ...
                'varAmp','obsAmp','period','varMax','obsMax','varMin','obsMin', ...
                'varPhase','obsPhase'};
        end

        function obj = fromSolverSnapshot(obs,snapshot)
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end

            obj = state();
            obj.obs = obs;
            obj.loadSolverSnapshot(snapshot);
        end

        function obj = fromSolverView(obs,solverView)
            solverView = fmam_state_ops.normalizeSolverView(solverView);
            derived = fmam_state_ops.buildDerivedView( ...
                obs,solverView,state.LphiConst,state.Lconst,state.countMax,state.errMax);
            obj = state.fromViews(obs,solverView,derived);
        end

        function obj = fromViews(obs,solverView,derived)
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end

            solverView = fmam_state_ops.normalizeSolverView(solverView);
            snapshot = fmam_state_ops.buildStateSnapshotFromViews(solverView,derived);
            obj = state.fromSolverSnapshot(obs,snapshot);
        end

        function [a,b] = primaryAmplitudeAndCenter(derived,PV)
            [a,b] = fmam_state_ops.primaryAmplitudeAndCenter(derived,PV);
        end
    end
end
