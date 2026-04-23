clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
dataDir = fullfile(scriptDir, 'data', 'fig3f');
saveFile = fullfile(dataDir, 'fig3f_data.mat');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(dataDir);

%% Temperature model and Fig. 3f targets
initialParams = [1.8746, 296.8523];
initialState = [1; 0];

startTarget = [80, 3];
orthMidTarget = [40, 3];
directMidTarget = [60, 4.5];
endTarget = [40, 6];

obs = [];
M = 50;
PV.name = 'var';
PV.idx = 2;

controlledIdx = [1, 2];
errBound = 1e-6;

continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);

orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 400], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

%% Build system and derivatives
run(fullfile(scriptDir, 'System_temp.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(initialParams));

%% Initial periodic orbit of the temperature model
initialOrbit = flexmod_forward_orbit(initialParams, initialState, struct('systemName', 'temp'));
if ~initialOrbit.success
    error('flexmod_refactor:Fig3fInitPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        initialOrbit.msg);
end

%% Start point
seedState = state(obs, initialParams, initialOrbit.orbit.t, initialOrbit.orbit.y, M, PV);
seedState.updatePeriod();
seedState.updateVar2();
seedSolverView = fmam_state_ops.solverViewFromState(seedState);

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = startTarget(1) / (2 * pi);
itemsPer(2).prop = 'varAmp';
itemsPer(2).idx = 2;
itemsPer(2).target = startTarget(2);

startTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
startTask.psiUpdateMode = true;
startTask.refreshPsiModeReferences();
startTask.needLog = false;

startTask.fit()
startTask.step()
startTask.errBound = 1e-12;
startTask.fit()

startDerivedView = startTask.exportDerivedView();
Parameters = reshape(startTask.exportSolverView().params, 1, []);
inferredTS = startDerivedView.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * startDerivedView.period, startDerivedView.t(end));
orbitOptions.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];

rhs1 = sys{1};
rhs2 = sys{2};
odeFunc = @(t, y, parameter) [ ...
    rhs1(reshape(y, 1, []), parameter); ...
    rhs2(reshape(y, 1, []), parameter)];
poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('flexmod_refactor:Fig3fStartPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

startPoint = struct();
startPoint.Parameters = Parameters;
startPoint.TS = {poResult.orbit_t, poResult.orbit_y};
startPoint.period = poFeatures.period;
startPoint.varAmp = reshape(poFeatures.state.amplitude, 1, []);
startPoint.varMax = reshape(poFeatures.state.max, 1, []);
startPoint.varMin = reshape(poFeatures.state.min, 1, []);
startPoint.logs = [];

%% Orthogonal period
seedState = state(obs, startPoint.Parameters, startPoint.TS{1}, startPoint.TS{2}, M, PV);
seedState.updatePeriod();
seedState.updateVar2();
seedSolverView = fmam_state_ops.solverViewFromState(seedState);

itemsPer(1).target = orthMidTarget(1) / (2 * pi);
itemsPer(2).target = orthMidTarget(2);

orthPeriodTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
orthPeriodTask.psiUpdateMode = true;
orthPeriodTask.refreshPsiModeReferences();
orthPeriodTask.needLog = true;

orthPeriodTask.fit()
orthPeriodTask.step()
orthPeriodTask.errBound = 1e-12;
orthPeriodTask.fit()

orthPeriodDerivedView = orthPeriodTask.exportDerivedView();
Parameters = reshape(orthPeriodTask.exportSolverView().params, 1, []);
inferredTS = orthPeriodDerivedView.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * orthPeriodDerivedView.period, orthPeriodDerivedView.t(end));
orbitOptions.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('flexmod_refactor:Fig3fOrthPeriodPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

orthPeriod = struct();
orthPeriod.Parameters = Parameters;
orthPeriod.TS = {poResult.orbit_t, poResult.orbit_y};
orthPeriod.period = poFeatures.period;
orthPeriod.varAmp = reshape(poFeatures.state.amplitude, 1, []);
orthPeriod.varMax = reshape(poFeatures.state.max, 1, []);
orthPeriod.varMin = reshape(poFeatures.state.min, 1, []);
orthPeriod.logs = orthPeriodTask.logs;

%% Orthogonal amplitude
seedState = state(obs, orthPeriod.Parameters, orthPeriod.TS{1}, orthPeriod.TS{2}, M, PV);
seedState.updatePeriod();
seedState.updateVar2();
seedSolverView = fmam_state_ops.solverViewFromState(seedState);

itemsPer(1).target = endTarget(1) / (2 * pi);
itemsPer(2).target = endTarget(2);

orthAmplitudeTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
orthAmplitudeTask.psiUpdateMode = true;
orthAmplitudeTask.refreshPsiModeReferences();
orthAmplitudeTask.needLog = true;

orthAmplitudeTask.fit()
orthAmplitudeTask.step()
orthAmplitudeTask.errBound = 1e-12;
orthAmplitudeTask.fit()

orthAmplitudeDerivedView = orthAmplitudeTask.exportDerivedView();
Parameters = reshape(orthAmplitudeTask.exportSolverView().params, 1, []);
inferredTS = orthAmplitudeDerivedView.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * orthAmplitudeDerivedView.period, orthAmplitudeDerivedView.t(end));
orbitOptions.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('flexmod_refactor:Fig3fOrthAmplitudePeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

orthAmplitude = struct();
orthAmplitude.Parameters = Parameters;
orthAmplitude.TS = {poResult.orbit_t, poResult.orbit_y};
orthAmplitude.period = poFeatures.period;
orthAmplitude.varAmp = reshape(poFeatures.state.amplitude, 1, []);
orthAmplitude.varMax = reshape(poFeatures.state.max, 1, []);
orthAmplitude.varMin = reshape(poFeatures.state.min, 1, []);
orthAmplitude.logs = orthAmplitudeTask.logs;

%% Direct middle
seedState = state(obs, startPoint.Parameters, startPoint.TS{1}, startPoint.TS{2}, M, PV);
seedState.updatePeriod();
seedState.updateVar2();
seedSolverView = fmam_state_ops.solverViewFromState(seedState);

itemsPer(1).target = directMidTarget(1) / (2 * pi);
itemsPer(2).target = directMidTarget(2);

directMidTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
directMidTask.psiUpdateMode = true;
directMidTask.refreshPsiModeReferences();
directMidTask.needLog = true;

directMidTask.fit()
directMidTask.step()
directMidTask.errBound = 1e-12;
directMidTask.fit()

directMidDerivedView = directMidTask.exportDerivedView();
Parameters = reshape(directMidTask.exportSolverView().params, 1, []);
inferredTS = directMidDerivedView.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * directMidDerivedView.period, directMidDerivedView.t(end));
orbitOptions.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('flexmod_refactor:Fig3fDirectMidPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

directMid = struct();
directMid.Parameters = Parameters;
directMid.TS = {poResult.orbit_t, poResult.orbit_y};
directMid.period = poFeatures.period;
directMid.varAmp = reshape(poFeatures.state.amplitude, 1, []);
directMid.varMax = reshape(poFeatures.state.max, 1, []);
directMid.varMin = reshape(poFeatures.state.min, 1, []);
directMid.logs = directMidTask.logs;

%% Direct end
seedState = state(obs, directMid.Parameters, directMid.TS{1}, directMid.TS{2}, M, PV);
seedState.updatePeriod();
seedState.updateVar2();
seedSolverView = fmam_state_ops.solverViewFromState(seedState);

itemsPer(1).target = endTarget(1) / (2 * pi);
itemsPer(2).target = endTarget(2);

directPathTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
directPathTask.psiUpdateMode = true;
directPathTask.refreshPsiModeReferences();
directPathTask.needLog = true;

directPathTask.fit()
directPathTask.step()
directPathTask.errBound = 1e-12;
directPathTask.fit()

directPathDerivedView = directPathTask.exportDerivedView();
Parameters = reshape(directPathTask.exportSolverView().params, 1, []);
inferredTS = directPathDerivedView.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * directPathDerivedView.period, directPathDerivedView.t(end));
orbitOptions.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('flexmod_refactor:Fig3fDirectEndPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

directPath = struct();
directPath.Parameters = Parameters;
directPath.TS = {poResult.orbit_t, poResult.orbit_y};
directPath.period = poFeatures.period;
directPath.varAmp = reshape(poFeatures.state.amplitude, 1, []);
directPath.varMax = reshape(poFeatures.state.max, 1, []);
directPath.varMin = reshape(poFeatures.state.min, 1, []);
directPath.logs = directPathTask.logs;

%% Save data
save(saveFile, ...
    'startPoint', ...
    'orthPeriod', ...
    'orthAmplitude', ...
    'directMid', ...
    'directPath');

fprintf('Saved data: %s\n', saveFile);
