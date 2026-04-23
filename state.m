classdef state < handle
    properties
        PV
        checkPsiNonnegative = true
    end

    properties (Access = private)
        timeMapWarningSignature = ''
        discretizationConfig = struct()
        extremaSearchConfig = struct()
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
        discretization
        extremaSearch

        t
        phi
        Psi
        TS_var
        TS_obs
    end
    methods
        function obj = state(obs,Params,t,TS_var,M,PV,varargin)
            obj.discretizationConfig = fmam_state_defaults.defaultDiscretization();
            obj.extremaSearchConfig = fmam_state_ops.defaultExtremaSearchSettings();
            if nargin == 0
                return
            end
            [obj.discretizationConfig,obj.extremaSearchConfig] = ...
                state.parseConstructorOptions(varargin{:});
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end
            t = t(:);
            Params = Params(:)';
            t = t - t(1);
            fmam_state_ops.validateTrajectoryInputs(obs,Params,t,TS_var,M,PV);

            obj.obs = obs;
            discretization = obj.discretization;
            extremaSearch = obj.extremaSearch;
            solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
                obs,Params,t,TS_var,M,PV,discretization,extremaSearch);
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
            discretization = obj.discretization;
            issue = fmam_state_ops.timeMapInvariantIssue( ...
                obj.Psi,obj.p_Psi,obj.q_Psi,discretization.reconstruction.phaseSampleCount);
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

        function prop = get.discretization(obj)
            if isempty(obj.discretizationConfig)
                obj.discretizationConfig = fmam_state_defaults.defaultDiscretization();
            end
            prop = obj.discretizationConfig;
        end

        function set.discretization(obj,discretization)
            current = obj.discretization;
            merged = state.mergeDiscretizationConfig(current,discretization);
            updated = fmam_state_defaults.normalizeDiscretization(merged);
            if isequal(updated,current)
                return
            end
            obj.discretizationConfig = updated;
            obj.refreshStateConfigurationCaches();
        end

        function prop = get.extremaSearch(obj)
            if isempty(obj.extremaSearchConfig)
                obj.extremaSearchConfig = fmam_state_ops.defaultExtremaSearchSettings();
            end
            prop = obj.extremaSearchConfig;
        end

        function set.extremaSearch(obj,extremaSearch)
            current = obj.extremaSearch;
            merged = state.mergeExtremaSearchConfig(current,extremaSearch);
            updated = fmam_state_ops.normalizeExtremaSearchSettings(merged);
            if isequal(updated,current)
                return
            end
            obj.extremaSearchConfig = updated;
            obj.refreshStateConfigurationCaches();
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
            L = obj.discretization.reconstruction.phaseSampleCount;
            prop = (0:L-1)'*2*pi/L;
        end

        function prop = get.t(obj)
            discretization = obj.discretization;
            issue = fmam_state_ops.timeMapInvariantIssue( ...
                obj.Psi,obj.p_Psi,obj.q_Psi,discretization.reconstruction.phaseSampleCount);
            obj.warnOnTimeMapIssue(issue);
            dt = fmam_state_ops.timeIncrementsFromCoefficients( ...
                obj.p_Psi,obj.q_Psi,discretization.reconstruction.phaseSampleCount);
            prop = [0; cumsum(dt(1:end-1))];
        end

        function updateVar2(obj)
            discretization = obj.discretization;
            extremaSearch = obj.extremaSearch;
            derived = fmam_state_ops.buildVariableDerivedState( ...
                obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi, ...
                discretization.reconstruction.phaseSampleCount, ...
                extremaSearch);

            obj.varPhiMax = derived.varPhiMax;
            obj.varPhiMin = derived.varPhiMin;
            obj.varAmp = derived.varAmp;
            obj.varMax = derived.varMax;
            obj.varMin = derived.varMin;
            obj.varPhase = derived.varPhase;
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

        function refreshDerivedState(obj)
            discretization = obj.discretization;
            extremaSearch = obj.extremaSearch;
            obj.assertTimeMapInvariant();
            varDerived = fmam_state_ops.buildVariableDerivedState( ...
                obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi, ...
                discretization.reconstruction.phaseSampleCount, ...
                extremaSearch);
            obsDerived = fmam_state_ops.buildObservableDerivedState( ...
                obj.obs,obj.p_var,obj.q_var,obj.p_Psi,obj.q_Psi,discretization,extremaSearch);

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
                obj.obs,solverView,discretization);

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
            [obj.a,obj.b] = fmam_state_ops.primaryAmplitudeAndCenter(derived,obj.PV);
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

        function refreshStateConfigurationCaches(obj)
            if state.hasLoadedSolverState(obj)
                obj.refreshDerivedState();
            end
        end

    end

    methods(Static)
        function obj = fromSolverSnapshot(obs,snapshot,discretization,extremaSearch)
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end

            obj = state();
            obj.obs = obs;
            if nargin < 3 || isempty(discretization)
                discretization = fmam_state_defaults.defaultDiscretization();
            end
            if nargin < 4 || isempty(extremaSearch)
                extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
            end
            obj.discretizationConfig = fmam_state_defaults.normalizeDiscretization(discretization);
            obj.extremaSearchConfig = fmam_state_ops.normalizeExtremaSearchSettings(extremaSearch);
            obj.loadSolverSnapshot(snapshot);
        end

        function obj = fromSolverView(obs,solverView,discretization,extremaSearch)
            solverView = fmam_state_ops.normalizeSolverView(solverView);
            if nargin < 3 || isempty(discretization)
                discretization = fmam_state_defaults.defaultDiscretization();
            end
            if nargin < 4 || isempty(extremaSearch)
                extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
            end
            derived = fmam_state_ops.buildDerivedView( ...
                obs,solverView,discretization);
            obj = state.fromViews(obs,solverView,derived,discretization,extremaSearch);
        end

        function obj = fromViews(obs,solverView,derived,discretization,extremaSearch)
            if isempty(obs)
                obs = {};
            elseif iscell(obs)
                obs = reshape(obs,1,[]);
            end

            solverView = fmam_state_ops.normalizeSolverView(solverView);
            snapshot = fmam_state_ops.buildStateSnapshotFromViews(solverView,derived);
            if nargin < 4 || isempty(discretization)
                discretization = fmam_state_defaults.defaultDiscretization();
            end
            if nargin < 5 || isempty(extremaSearch)
                extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
            end
            obj = state.fromSolverSnapshot(obs,snapshot,discretization,extremaSearch);
        end
    end

    methods(Static, Access = private)
        function [discretization,extremaSearch] = parseConstructorOptions(varargin)
            discretization = fmam_state_defaults.defaultDiscretization();
            extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
            if isempty(varargin)
                return
            end

            if numel(varargin) <= 2 && all(cellfun(@isstruct,varargin))
                discretization = fmam_state_defaults.normalizeDiscretization(varargin{1});
                if numel(varargin) == 2
                    extremaSearch = fmam_state_ops.normalizeExtremaSearchSettings(varargin{2});
                end
                return
            end

            if mod(numel(varargin),2) ~= 0
                error('state:InvalidConstructorOptions', ...
                    'Optional state constructor arguments must be structs or name-value pairs.')
            end
            for i = 1:2:numel(varargin)
                name = varargin{i};
                value = varargin{i + 1};
                if ~ischar(name) && ~(isstring(name) && isscalar(name))
                    error('state:InvalidConstructorOptions', ...
                        'Option names must be character vectors or scalar strings.')
                end
                switch lower(char(name))
                    case 'discretization'
                        discretization = fmam_state_defaults.normalizeDiscretization(value);
                    case 'extremasearch'
                        extremaSearch = fmam_state_ops.normalizeExtremaSearchSettings(value);
                    otherwise
                        error('state:InvalidConstructorOptions', ...
                            'Unsupported state constructor option ''%s''.',char(name))
                end
            end
        end

        function discretization = mergeDiscretizationConfig(current,update)
            if nargin < 1 || isempty(current)
                current = fmam_state_defaults.defaultDiscretization();
            end
            if nargin < 2 || isempty(update)
                discretization = current;
                return
            end
            if ~isstruct(update)
                error('state:InvalidDiscretizationUpdate', ...
                    'discretization updates must be provided as a struct.');
            end

            discretization = current;
            fieldNames = fieldnames(update);
            for i = 1:numel(fieldNames)
                fieldName = fieldNames{i};
                if strcmp(fieldName,'reconstruction')
                    discretization.reconstruction = fmam_state_defaults.normalizeReconstruction( ...
                        update.reconstruction, current.reconstruction);
                else
                    discretization.(fieldName) = update.(fieldName);
                end
            end
        end

        function extremaSearch = mergeExtremaSearchConfig(current,update)
            if nargin < 1 || isempty(current)
                current = fmam_state_ops.defaultExtremaSearchSettings();
            end
            if nargin < 2 || isempty(update)
                extremaSearch = current;
                return
            end
            if ~isstruct(update)
                error('state:InvalidExtremaSearchUpdate', ...
                    'extremaSearch updates must be provided as a struct.');
            end

            extremaSearch = current;
            fieldNames = fieldnames(update);
            for i = 1:numel(fieldNames)
                extremaSearch.(fieldNames{i}) = update.(fieldNames{i});
            end
        end

        function tf = hasLoadedSolverState(obj)
            tf = ~isempty(obj.p_Psi) && ~isempty(obj.p_var) && ~isempty(obj.q_var) && ...
                ~isempty(obj.params) && ~isempty(obj.PV);
        end

        function fieldNames = solverStateFields()
            fieldNames = { ...
                'PV','params','p_Psi','q_Psi','p_var','q_var', ...
                'p_var_origin','q_var_origin','p_obs','q_obs', ...
                'a','b','varPhiMax','obsPhiMax','varPhiMin','obsPhiMin', ...
                'varAmp','obsAmp','period','varMax','obsMax','varMin','obsMin', ...
                'varPhase','obsPhase'};
        end
    end
end
