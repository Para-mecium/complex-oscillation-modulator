clear
clc
%%
scriptDir = fileparts(mfilename('fullpath'));
circuitDir = fileparts(scriptDir);
fmamDir = fileparts(circuitDir);

multiplier = 2;
needPath = true;
%%
rawData = load(fullfile(circuitDir, 'learnedData_ODE.mat'));

switch multiplier
    case 1
        targetTag = '1p0x';
    case 1.5
        targetTag = '1p5x';
    case 2
        targetTag = '2p0x';
end
targetfile = fullfile(scriptDir, ['period_target_' targetTag '.mat']);

M = 75;
PV = struct(...
    'name', 'var',...
    'idx', 1);
params = reshape(rawData.Parameters, 1, []);

t = rawData.TS{1};
TS_var = rawData.TS{2};
obs = [];
%%
startDir = pwd;
cleanupObj = onCleanup(@() cd(startDir));

cd(circuitDir);
System

cd(fmamDir);
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

stateObj = state(obs, params, t, TS_var, M, PV);
stateObj.updatePeriod();
stateObj.updateVar2();

solverView = fmam_state_ops.solverViewFromState(stateObj);

basePeriod = stateObj.period;
baseVarAmp1 = stateObj.varAmp(1);
periodTarget = basePeriod * multiplier;

items_per = struct([]);
items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = periodTarget / (2 * pi);

items_per(2).prop = 'varAmp';
items_per(2).idx = 1;
items_per(2).target = baseVarAmp1;

items_controlled = [1, 4];

errBound = 1e-6;
continuationOptions = struct('initialLambdaStep', 0.01);

task = FMAM_ODE(sys, obs, solverView, items_per, items_controlled, ...
    [], errBound, 'derivatives', derivatives, ...
    'continuationOptions', continuationOptions);
task.isPsiUpdated = true;

task.fit();
task.step();
task.errBound = 1e-12;
task.fit();

stateView = task.exportSolverView();
stateDerived = task.exportDerivedView();

Parameters = reshape(stateView.params, 1, []);

%%
addpath(fullfile(fmamDir, 'PO_extract'));
odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter)];
y0 = stateDerived.TS_var(1, :).';
searchWindow = max(10 * stateDerived.period, stateDerived.t(end));
orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, searchWindow], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('Modulation_bi_target:PeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction failed for %.1fx (%s).', ...
        multiplier, poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
features = evaluate_orbit_features(orbitForFeatures, [], [], struct());

nSeries = size(features.state.series, 1);
TS = {poResult.orbit_t(1:nSeries), features.state.series};

period = features.period;
varAmp = reshape(features.state.amplitude, 1, []);
varMax = reshape(features.state.max, 1, []);
varMin = reshape(features.state.min, 1, []);
periodMultiplier = multiplier;

save(targetfile, ...
    'Parameters', 'TS', 'period', 'varAmp', 'varMax', 'varMin', ...
    'periodMultiplier', 'periodTarget');

fprintf('Saved %s\n', targetfile);

%%
if needPath
    lambdaStepCap = 5e-3;
    pathFile = fullfile(scriptDir, 'params_modulation_path.mat');
    continuationOptionsPath = struct( ...
        'initialLambdaStep', lambdaStepCap, ...
        'maxLambdaStep', lambdaStepCap);
    taskPath = FMAM_ODE(sys, obs, solverView, items_per, items_controlled, ...
        [], errBound, 'derivatives', derivatives, ...
        'continuationOptions', continuationOptionsPath);
    taskPath.isPsiUpdated = true;
    taskPath.needLog = true;
    taskPath.fit();
    params_start = reshape(taskPath.exportSolverView().params, 1, []);
    taskPath.step();
    taskPath.errBound = 1e-12;
    taskPath.fit();

    stateViewPath = taskPath.exportSolverView();
    params_end = reshape(double(stateViewPath.params), 1, []);
    solution = taskPath.logs;

    curve_params = zeros(numel(solution), numel(solution(1).params));
    for i = 1:numel(solution)
        curve_params(i, :) = reshape(double(solution(i).params), 1, []);
    end
    pathData = struct();
    pathData.curve_params = curve_params;
    pathData.params_start = params_start;
    pathData.params_end = params_end;
    save(pathFile, '-struct', 'pathData');
    fprintf('Saved %s\n', pathFile);
end
