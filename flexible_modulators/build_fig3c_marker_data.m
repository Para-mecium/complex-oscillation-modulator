clear
clc

% Three Fig. 3c markers: (80, 2.0), (80, 2.5), (80, 3.0)
targetPeriod = 80; 
targetAmplitude = 3.0;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'data', 'common', 'initData.mat');
markerDataDir = fullfile(scriptDir, 'data', 'fig3c', 'markers');

amplitudeTag = strrep(sprintf('%.1f', targetAmplitude), '.', 'p');
saveFile = fullfile(markerDataDir, ...
    sprintf('fig3c_marker_T%03d_A%s.mat', round(targetPeriod), amplitudeTag));

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(markerDataDir);
%% Single marker settings for Fig. 3c
initData = load(initDataFile);
startParams = initData.Parameters;
startT = initData.TS{1};
startTSVar = initData.TS{2};

obs = [];
M = 50;
PV.name = 'var';
PV.idx = 2;

controlledIdx = [1, 2];
errBound = 1e-6;

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = targetPeriod / (2 * pi);

itemsPer(2).prop = 'varAmp';
itemsPer(2).idx = 2;
itemsPer(2).target = targetAmplitude;

continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);

%% Build system and derivatives
run(fullfile(scriptDir, 'System_base.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(startParams));

%% Build FMAM state from initData
seedState = state(obs, startParams, startT, startTSVar, M, PV);
seedState.updatePeriod();
seedState.updateVar2();
seedSolverView = fmam_state_ops.solverViewFromState(seedState);

%% Compute the marker point
markerTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
markerTask.isPsiUpdated = true;
markerTask.needLog = false;

markerTask.fit()
markerTask.step()
markerTask.errBound = 1e-12;
markerTask.fit()

markerDerivedView = markerTask.exportDerivedView();
Parameters = reshape(markerTask.exportSolverView().params, 1, []);

%% Forward integration and feature extraction at the marker point
inferredTS = markerDerivedView.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * markerDerivedView.period, markerDerivedView.t(end));
orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 400], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);
orbitOptions.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];

rhs1 = sys{1};
rhs2 = sys{2};
odeFunc = @(t, y, parameter) [ ...
    rhs1(reshape(y, 1, []), parameter); ...
    rhs2(reshape(y, 1, []), parameter)];
poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('flexmod_refactor:Fig3cMarkerPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

TS = {poResult.orbit_t, poResult.orbit_y};
period = poFeatures.period;
varAmp = reshape(poFeatures.state.amplitude, 1, []);
varMax = reshape(poFeatures.state.max, 1, []);
varMin = reshape(poFeatures.state.min, 1, []);

save(saveFile, ...
    'Parameters', ...
    'TS', ...
    'period', ...
    'varAmp', ...
    'varMax', ...
    'varMin');

fprintf('Saved marker data: %s\n', saveFile);
