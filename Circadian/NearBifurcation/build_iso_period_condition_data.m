clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
circadianDir = fileparts(scriptDir);
repoDir = fileparts(circadianDir);
initDataFile = fullfile(circadianDir, 'data', 'common', 'initData.mat');
curveDataDir = fullfile(scriptDir, 'data', 'iso_period', 'curves');

addpath(repoDir);
addpath(circadianDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(curveDataDir);

%% Iso-period settings
% targetPeriods = [23.5, 24.0, 24.5, 25.0];
% seedAt = [0.070, 0.060, 0.050, 0.040];
targetPeriods = [25.0];
seedAt = [0.040];
downAt = 0;
upAt = 0.09;
curveStepCap = 1e-2;

obs = {@(variable) variable(:, 2) + variable(:, 3)};
M = 50;
PV.name = 'obs';
PV.idx = 1;

controlledIdx = [1, 2];
errBound = 1e-6;
conditionStopRcond = 1e-12;

%% Build system and derivatives
initData = load(initDataFile);
baseInitialParams = initData.Parameters;
baseT = initData.TS{1};
baseTSVar = initData.TS{2};
baseY0 = initData.TS{2}(1, :).';

run(fullfile(circadianDir, 'System.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(baseInitialParams));

%% Build base FMAM state from initData
baseState = state(obs, baseInitialParams, baseT, baseTSVar, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);

%% Compute one seed and two A_T branches for each target period
for i = 1:numel(targetPeriods)
    targetPeriod = targetPeriods(i);
    targetSeedAt = seedAt(i);
    outputFile = fullfile(curveDataDir, ...
        sprintf('iso_period_T%s_condition.mat', period_tag(targetPeriod)));

    seedTask = solve_seed_task(baseSolverView, targetPeriod, targetSeedAt, controlledIdx, errBound, derivatives, sys, conditionStopRcond);
    seed = export_seed_data(seedTask, baseY0);
    seedStatus = seedTask.continuationStatus;

    downBranchTask = solve_branch_task( ...
        seedTask.exportSolverView(), targetPeriod, downAt, curveStepCap, controlledIdx, errBound, derivatives, sys, conditionStopRcond);
    upBranchTask = solve_branch_task( ...
        seedTask.exportSolverView(), targetPeriod, upAt, curveStepCap, controlledIdx, errBound, derivatives, sys, conditionStopRcond);

    downBranch = downBranchTask.logs;
    upBranch = upBranchTask.logs;
    downStatus = downBranchTask.continuationStatus;
    upStatus = upBranchTask.continuationStatus;

    save(outputFile, ...
        'targetPeriod', ...
        'targetSeedAt', ...
        'downAt', ...
        'upAt', ...
        'conditionStopRcond', ...
        'seed', ...
        'seedStatus', ...
        'downBranch', ...
        'upBranch', ...
        'downStatus', ...
        'upStatus');

    fprintf('Saved data: %s\n', outputFile);
end

%%
function task = solve_seed_task(baseSolverView, targetPeriod, targetAt, controlledIdx, errBound, derivatives, sys, conditionStopRcond)
continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', conditionStopRcond);

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = targetPeriod / (2 * pi);

itemsPer(2).prop = 'params';
itemsPer(2).idx = 2;
itemsPer(2).target = targetAt;

task = FMAM_ODE(sys, obs_spec(), baseSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
task.psiUpdateMode = true;
task.refreshPsiModeReferences();
task.needLog = true;

task.fit()
task.step()
task.errBound = 1e-12;
task.fit()
end

%%
function task = solve_branch_task(seedSolverView, targetPeriod, targetAt, curveStepCap, controlledIdx, errBound, derivatives, sys, conditionStopRcond)
continuationOptions = struct( ...
    'initialLambdaStep', min(curveStepCap, 0.05), ...
    'maxLambdaStep', curveStepCap, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', conditionStopRcond);

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = targetPeriod / (2 * pi);

itemsPer(2).prop = 'params';
itemsPer(2).idx = 2;
itemsPer(2).target = targetAt;

task = FMAM_ODE(sys, obs_spec(), seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
task.psiUpdateMode = true;
task.refreshPsiModeReferences();
task.needLog = true;

task.fit()
task.step()
end

function seedData = export_seed_data(task, y0)
finalParams = reshape(task.exportSolverView().params, 1, []);
orbitResult = circadian_forward_orbit(finalParams, y0, struct());
if ~orbitResult.success
    error('circadian_refactor:IsoPeriodSeedOrbitGenerationFailed', ...
        'Periodic-orbit extraction failed at seed point.');
end

seedData = struct();
seedData.Parameters = finalParams;
seedData.TS = {orbitResult.orbit.t, orbitResult.orbit.y};
seedData.period = orbitResult.features.period;
seedData.obsAmp = reshape(orbitResult.features.observable.amplitude, 1, []);
seedData.obsMax = reshape(orbitResult.features.observable.max, 1, []);
seedData.obsMin = reshape(orbitResult.features.observable.min, 1, []);
seedData.directConditionEstimate = task.logs(end).directConditionEstimate;
seedData.conditionNumber = 1 / seedData.directConditionEstimate;
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end

function obs = obs_spec()
obs = {@(variable) variable(:, 2) + variable(:, 3)};
end
