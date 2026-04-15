function result = run_flexmod_fmam_task(cfg)
if nargin < 1
    cfg = struct();
end

flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);
model = flexmod.build_model(cfg);

if isfield(cfg, 'initialSolverView') && ~isempty(cfg.initialSolverView)
    currentView = normalize_solver_view_input(cfg.initialSolverView);
    solverInput = currentView;
    currentDerived = resolve_initial_derived_view(cfg, model, currentView);
    startParams = reshape(currentView.params, 1, []);
    seedOrbit = flexmod.orbit_from_state(currentDerived, startParams);
elseif isfield(cfg, 'initialStat') && isa(cfg.initialStat, 'state')
    stat = cfg.initialStat;
    startParams = reshape(stat.params, 1, []);
    seedOrbit = flexmod.orbit_from_state(stat, startParams);
    currentView = fmam_state_ops.solverViewFromState(stat);
    solverInput = currentView;
    currentDerived = fmam_state_ops.derivedViewFromState(stat);
else
    if isempty(cfg.startParams)
        startParams = model.defaultParams;
    else
        startParams = reshape(cfg.startParams, 1, []);
    end

    seedOrbit = flexmod.find_orbit(model, startParams, cfg, cfg.initialState);
    if ~seedOrbit.success
        error('flexmod:SeedOrbitFailed', 'Failed to locate a seed periodic orbit: %s', seedOrbit.message);
    end

    stat = state(model.observables, startParams, seedOrbit.t, seedOrbit.y, cfg.truncationOrder, model.pv);
    stat.updatePeriod();
    stat.updateVar2();
    currentView = fmam_state_ops.solverViewFromState(stat);
    solverInput = currentView;
    currentDerived = fmam_state_ops.derivedViewFromState(stat);
end

[itemsPerturb, itemsControlled, maxStep] = build_goal_system(cfg, model, currentView, currentDerived);
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
% task.verbose = logical(cfg.fmam.verbose);

task.fit();
task.step();

finalView = task.exportSolverView();
finalDerived = task.exportDerivedView();
finalParams = reshape(finalView.params, 1, []);
finalOrbit = flexmod.refine_orbit_from_state(model, finalParams, cfg, finalDerived);
path = flexmod.reconstruct_path_data(model, startParams, seedOrbit, task.logs, task.continuationStatus, cfg);

result = struct();
result.cfg = cfg;
result.model = model;
result.params = finalParams;
result.orbit = finalOrbit;
result.measures = struct( ...
    'period', finalOrbit.period, ...
    'amplitude', finalOrbit.amplitude, ...
    'y_max', finalOrbit.yMax, ...
    'y_min', finalOrbit.yMin);
result.path = path;
result.logs = task.logs;
result.fmamStatus = task.continuationStatus;
result.fitResult = task;
result.solverView = finalView;
result.derivedView = finalDerived;
end

function [itemsPerturb, itemsControlled, maxStep] = build_goal_system(cfg, model, solverView, derivedView)
goalOrder = cfg.goalOrder;
itemsPerturb = repmat(struct('prop', '', 'idx', 1, 'target', NaN), 1, numel(goalOrder));
itemsControlled = zeros(1, numel(goalOrder));
delta = zeros(1, numel(goalOrder));

for i = 1:numel(goalOrder)
    name = goalOrder{i};
    switch lower(name)
        case 'period'
            itemsPerturb(i).prop = 'p_Psi';
            itemsPerturb(i).idx = 1;
            itemsPerturb(i).target = cfg.goals.period / (2 * pi);
            delta(i) = itemsPerturb(i).target - solverView.p_Psi(1);
        case 'amplitude'
            itemsPerturb(i).prop = 'varAmp';
            itemsPerturb(i).idx = 2;
            itemsPerturb(i).target = cfg.goals.amplitude;
            delta(i) = itemsPerturb(i).target - derivedView.varAmp(2);
        otherwise
            idx = find(strcmp(model.paramNames, name), 1, 'first');
            itemsPerturb(i).prop = 'params';
            itemsPerturb(i).idx = idx;
            itemsPerturb(i).target = cfg.goals.(name);
            delta(i) = itemsPerturb(i).target - solverView.params(idx);
    end

    controlName = cfg.controlledParams{i};
    itemsControlled(i) = find(strcmp(model.paramNames, controlName), 1, 'first');
end

maxStep = max(abs(delta) / cfg.fmam.targetSteps, cfg.fmam.minGoalStep);
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
