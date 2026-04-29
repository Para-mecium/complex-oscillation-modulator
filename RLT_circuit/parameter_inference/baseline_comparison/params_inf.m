clear
clc

%% Proposed FMAM parameter inference
needPath = true;
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');
addpath(scriptDir, '-begin');

outputDir = fullfile(scriptDir, 'results', 'proposed_method');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

targetDataFile = fullfile(parameterInferenceDir, 'initData_circuit.mat');
baseDataFile = fullfile(parameterInferenceDir, 'initData_ODE.mat');

targetData = load(targetDataFile);
initData = load(baseDataFile);

run(fullfile(circuitDir, 'System.m'));
obs = [];

params = initData.Parameters;
t = initData.TS{1};
TS_var = initData.TS{2};
M = 50;
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

PV.name = 'var';
PV.idx = 1;
State = state(obs, params, t, TS_var, M, PV);
StateView = fmam_state_ops.solverViewFromState(State);

items_per = struct;
items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = targetData.period / (2 * pi);

items_per(2).prop = 'varAmp';
items_per(2).idx = 1;
items_per(2).target = targetData.varAmp(1);

items_per(3).prop = 'varAmp';
items_per(3).idx = 2;
items_per(3).target = targetData.varAmp(2);

items_per(4).prop = 'varAmp';
items_per(4).idx = 3;
items_per(4).target = targetData.varAmp(3);

items_controlled = [1 4 5 6];
errBound = 1e-6;
continuationOptions = struct('initialLambdaStep', 0.01, 'predictorMode', 'constant');
Modtask = FMAM_ODE(sys, obs, StateView, items_per, items_controlled, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
Modtask.psiUpdateMode = true;
Modtask.refreshPsiModeReferences();
Modtask.needLog = needPath;

runTimer = tic;
Modtask.fit()
Modtask.step()
Modtask.errBound = 1e-12;
Modtask.fit()
Modtask.step()
elapsedTime = toc(runTimer);
fprintf('Proposed method computing time: %.6f seconds\n', elapsedTime);

StateView = Modtask.exportSolverView();
StateDerived = Modtask.exportDerivedView();

%% Extract ODE periodic orbit at inferred parameters
Parameters = reshape(double(StateView.params), 1, []);
inferredTS = StateDerived.TS_var;
y0 = inferredTS(1, :).';
odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter)];
searchWindow = max(10 * StateDerived.period, StateDerived.t(end));
poExtractDir = fullfile(repoDir, 'PO_extract');
addpath(poExtractDir, '-begin');

poOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, searchWindow], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, poOptions);

if ~poResult.has_orbit
    error('params_inf:PeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbit = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbit, [], [], struct());

%% Save proposed-method result in baseline-comparison format
userConfig = struct();
userConfig.targetDataFile = targetDataFile;
userConfig.baseDataFile = baseDataFile;
userConfig.forwardOptions.poOptions.solver_tol = poOptions.solver_tol;
config = build_config(userConfig);

bestLoss = loss_function(orbit, config.targetOrbit, config.lossOptions);
activeParams = Parameters(config.activeIndex);

result = struct();
result.methodName = 'proposed_method';
result.bestLoss = bestLoss;
result.bestActiveParams = activeParams;
result.bestParameters = Parameters;
result.runtime = elapsedTime;
result.totalBestSoFarLosses = bestLoss;
result.totalTimeSoFar = elapsedTime;
result.bestForwardResult = struct( ...
    'success', true, ...
    'msg', "", ...
    'Parameters', Parameters, ...
    'initialState', y0, ...
    'orbit', orbit, ...
    'features', poFeatures);
result.targetDataFile = targetDataFile;
result.baseDataFile = baseDataFile;
result.activeIndex = config.activeIndex;
result.activeNames = config.activeNames;
result.activeLowerBound = config.activeLowerBound;
result.activeUpperBound = config.activeUpperBound;
result.controlledIndex = items_controlled;
result.items = items_per;
result.errBound = Modtask.errBound;
result.needPath = needPath;

if needPath
    solutionPath = Modtask.logs;
    curveParams = zeros(numel(solutionPath), numel(Parameters));
    for i = 1:numel(solutionPath)
        curveParams(i, :) = reshape(double(solutionPath(i).params), 1, []);
    end
    result.curveParams = curveParams;
    result.startParameters = reshape(double(State.params), 1, []);
end

configSummary = config_summary(config);
resultFile = fullfile(outputDir, 'proposed_method_result.mat');
save(resultFile, 'result', 'configSummary');
fprintf('Saved proposed-method baseline-comparison result: %s\n', resultFile);
