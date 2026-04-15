function result = run_circadian_fmam_task(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
model = circadian.build_model(cfg);

if isfield(cfg, 'initialSolverView') && ~isempty(cfg.initialSolverView)
    currentView = normalize_solver_view_input(cfg.initialSolverView);
    solverInput = currentView;
    currentDerived = resolve_initial_derived_view(cfg, model, currentView);
    startParams = reshape(currentView.params, 1, []);
    seedOrbit = circadian.orbit_from_state(currentDerived, startParams);
elseif isfield(cfg, 'initialStat') && isa(cfg.initialStat, 'state')
    stat = cfg.initialStat;
    startParams = reshape(stat.params, 1, []);
    seedOrbit = circadian.orbit_from_state(stat, startParams);
    currentView = fmam_state_ops.solverViewFromState(stat);
    solverInput = currentView;
    currentDerived = fmam_state_ops.derivedViewFromState(stat);
else
    if isempty(cfg.startParams)
        startParams = model.defaultParams;
    else
        startParams = reshape(cfg.startParams, 1, []);
    end

    seedOrbit = circadian.find_orbit(model, startParams, cfg, cfg.initialState);
    if ~seedOrbit.success
        error('circadian:SeedOrbitFailed', ...
            'Failed to locate a seed periodic orbit: %s', seedOrbit.message);
    end

    stat = state(model.observables, startParams, seedOrbit.t, seedOrbit.y, cfg.truncationOrder, model.pv);
    stat.updatePeriod();
    stat.updateVar2();
    currentView = fmam_state_ops.solverViewFromState(stat);
    solverInput = currentView;
    currentDerived = fmam_state_ops.derivedViewFromState(stat);
end

[itemsPerturb, itemsControlled, maxStep] = circadian.build_task_inputs( ...
    cfg, model, currentView, currentDerived, cfg);
continuationOptions = struct( ...
    'predictorMode', cfg.fmam.predictorMode, ...
    'initialLambdaStep', min(cfg.fmam.dlambdaCap, 0.05), ...
    'maxLambdaStep', cfg.fmam.dlambdaCap, ...
    'conditionStopEnabled', cfg.fmam.conditionStopEnabled, ...
    'conditionStopRcond', cfg.fmam.conditionStopRcond);

task = FMAM_ODE(model.system, model.observables, solverInput, itemsPerturb, itemsControlled, ...
    maxStep, cfg.fmam.errBound, 'derivatives', model.derivatives, ...
    'continuationOptions', continuationOptions);
task.isPsiUpdated = true;
task.needLog = true;
task.verbose = logical(cfg.fmam.verbose);

task.fit();
task.step();

finalView = task.exportSolverView();
finalDerived = task.exportDerivedView();
finalParams = reshape(finalView.params, 1, []);
finalOrbit = circadian.orbit_from_state(finalDerived, finalParams);
path = circadian.reconstruct_path_data(model, startParams, seedOrbit, task.logs, task.continuationStatus, cfg);

result = struct();
result.cfg = cfg;
result.model = model;
result.params = finalParams;
result.orbit = finalOrbit;
result.measures = struct( ...
    'period', finalDerived.period, ...
    'amplitude', finalDerived.obsAmp(1), ...
    'maximum', finalDerived.obsMax(1), ...
    'minimum', finalDerived.obsMin(1));
result.path = path;
result.logs = task.logs;
result.fmamStatus = task.continuationStatus;
result.fitResult = task;
result.solverView = finalView;
result.derivedView = finalDerived;
end

function derivedView = resolve_initial_derived_view(cfg, model, solverView)
if isfield(cfg, 'initialDerivedView') && ~isempty(cfg.initialDerivedView)
    derivedView = cfg.initialDerivedView;
else
    discretization = fmam_state_defaults.defaultDiscretization();
    derivedView = fmam_state_ops.buildDerivedView( ...
        model.observables, solverView, discretization);
end
end

function solverView = normalize_solver_view_input(view)
if isstruct(view) && isfield(view, 'solverView')
    solverView = view.solverView;
else
    solverView = view;
end
solverView = fmam_state_ops.normalizeSolverView(solverView);
end
