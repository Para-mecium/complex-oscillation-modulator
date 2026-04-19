clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'data', 'common', 'initData.mat');
curveDataDir = fullfile(scriptDir, 'data', 'fig5b', 'curves');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(curveDataDir);

%% Fig. 5b curve settings
targetPeriods = [23.5, 24.0, 24.5, 25.0];
seedKd = 1e-4 * [1, 0.5, 0.5, 0.25];
lowerSeedAt = [0.04, 0.035, 0.030, 0.02];
upperSeedAt = [0.070, 0.060, 0.050, 0.040];

leftKd = 0;
rightKd = 1.5e-4;
curveStepCap = 5e-3;

obs = {@(variable) variable(:, 2) + variable(:, 3)};
M = 50;
PV.name = 'obs';
PV.idx = 1;

controlledIdx = [1, 2];
errBound = 1e-6;

%% Build system and derivatives
initData = load(initDataFile);
baseInitialParams = initData.Parameters;
baseT = initData.TS{1};
baseTSVar = initData.TS{2};
baseY0 = initData.TS{2}(1, :).';

run(fullfile(scriptDir, 'System.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(baseInitialParams));

%% Build base FMAM state from initData
baseState = state(obs, baseInitialParams, baseT, baseTSVar, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);

%% Compute two seeds and four branches for each target period
for i = 1:numel(targetPeriods)
    targetPeriod = targetPeriods(i);
    lowerAt = lowerSeedAt(i);
    upperAt = upperSeedAt(i);
    outputFile = fullfile(curveDataDir, ...
        sprintf('fig5b_iso_period_curve_T%s.mat', period_tag(targetPeriod)));

    lowerSeedTask = solve_seed_task(baseSolverView, targetPeriod, lowerAt, controlledIdx, errBound, derivatives, sys);
    upperSeedTask = solve_seed_task(baseSolverView, targetPeriod, upperAt, controlledIdx, errBound, derivatives, sys);

    lowerSeed = export_seed_data(lowerSeedTask, baseY0);
    upperSeed = export_seed_data(upperSeedTask, baseY0);

    lowerLeftBranchTask = solve_branch_task( ...
        lowerSeedTask.exportSolverView(), targetPeriod, leftKd, curveStepCap, controlledIdx, errBound, derivatives, sys);
    lowerRightBranchTask = solve_branch_task( ...
        lowerSeedTask.exportSolverView(), targetPeriod, rightKd, curveStepCap, controlledIdx, errBound, derivatives, sys);
    upperLeftBranchTask = solve_branch_task( ...
        upperSeedTask.exportSolverView(), targetPeriod, leftKd, curveStepCap, controlledIdx, errBound, derivatives, sys);
    upperRightBranchTask = solve_branch_task( ...
        upperSeedTask.exportSolverView(), targetPeriod, rightKd, curveStepCap, controlledIdx, errBound, derivatives, sys);

    lowerLeftBranch = lowerLeftBranchTask.logs;
    lowerRightBranch = lowerRightBranchTask.logs;
    upperLeftBranch = upperLeftBranchTask.logs;
    upperRightBranch = upperRightBranchTask.logs;

    save(outputFile, ...
        'targetPeriod', ...
        'lowerSeed', ...
        'upperSeed', ...
        'lowerLeftBranch', ...
        'lowerRightBranch', ...
        'upperLeftBranch', ...
        'upperRightBranch');

    fprintf('Saved data: %s\n', outputFile);
end

%%
function task = solve_seed_task(baseSolverView, targetPeriod, targetAt, controlledIdx, errBound, derivatives, sys)
continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = targetPeriod / (2 * pi);

itemsPer(2).prop = 'params';
itemsPer(2).idx = 2;
itemsPer(2).target = targetAt;

task = FMAM_ODE(sys, obs_spec(), baseSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
task.isPsiUpdated = true;
task.needLog = false;

task.fit()
task.step()
task.errBound = 1e-12;
task.fit()
end

%%
function task = solve_branch_task(seedSolverView, targetPeriod, targetKd, curveStepCap, controlledIdx, errBound, derivatives, sys)
continuationOptions = struct( ...
    'initialLambdaStep', min(curveStepCap, 0.05), ...
    'maxLambdaStep', curveStepCap, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);

itemsPer = struct([]);
itemsPer(1).prop = 'p_Psi';
itemsPer(1).idx = 1;
itemsPer(1).target = targetPeriod / (2 * pi);

itemsPer(2).prop = 'params';
itemsPer(2).idx = 1;
itemsPer(2).target = targetKd;

task = FMAM_ODE(sys, obs_spec(), seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
task.isPsiUpdated = true;
task.needLog = true;

task.fit()
task.step()
end

function seedData = export_seed_data(task, y0)
finalParams = reshape(task.exportSolverView().params, 1, []);
orbitResult = circadian_forward_orbit(finalParams, y0, struct());
if ~orbitResult.success
    error('circadian_refactor:Fig5bSeedOrbitGenerationFailed', ...
        'Periodic-orbit extraction failed at seed point.');
end

seedData = struct();
seedData.Parameters = finalParams;
seedData.TS = {orbitResult.orbit.t, orbitResult.orbit.y};
seedData.period = orbitResult.features.period;
seedData.obsAmp = reshape(orbitResult.features.observable.amplitude, 1, []);
seedData.obsMax = reshape(orbitResult.features.observable.max, 1, []);
seedData.obsMin = reshape(orbitResult.features.observable.min, 1, []);
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end

function obs = obs_spec()
obs = {@(variable) variable(:, 2) + variable(:, 3)};
end
