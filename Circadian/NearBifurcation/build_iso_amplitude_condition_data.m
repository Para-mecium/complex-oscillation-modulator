clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
circadianDir = fileparts(scriptDir);
repoDir = fileparts(circadianDir);
initDataFile = fullfile(circadianDir, 'data', 'common', 'initData.mat');
curveDataDir = fullfile(scriptDir, 'data', 'iso_amplitude', 'curves');

addpath(repoDir);
addpath(circadianDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(curveDataDir);

%% Iso-amplitude settings
targetAmplitudes = [0.005, 0.015, 0.025, 0.035, 0.0385];
centerPeriod = 24.0;

leftKd = 0;
rightKd = 1.5e-4;

obs = {@(variable) variable(:, 2) + variable(:, 3)};
M = 50;
PV.name = 'obs';
PV.idx = 1;

controlledIdx = [1, 2];
errBound = 1e-6;
curveStepCap = 1e-2;
conditionStopRcond = 1e-10;

%% Build system and derivatives
run(fullfile(circadianDir, 'System.m'));
derivatives = build_symbolic_derivatives(sys, obs, 2);

%% Load init orbit and build base state
initData = load(initDataFile);
baseInitialParams = initData.Parameters;
baseT = initData.TS{1};
baseTSVar = initData.TS{2};
baseY0 = initData.TS{2}(1, :).';

baseState = state(obs, baseInitialParams, baseT, baseTSVar, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);

%% Compute one seed and two branches for each target amplitude
for i = 1:numel(targetAmplitudes)
    targetAmplitude = targetAmplitudes(i);

    seedTask = solve_seed_task(baseSolverView, targetAmplitude, centerPeriod, controlledIdx, errBound, derivatives, sys, conditionStopRcond);
    seed = export_seed_data(seedTask, baseY0);
    seedStatus = seedTask.continuationStatus;

    leftBranchTask = solve_branch_task( ...
        seedTask.exportSolverView(), targetAmplitude, leftKd, curveStepCap, controlledIdx, errBound, derivatives, sys, conditionStopRcond);
    rightBranchTask = solve_branch_task( ...
        seedTask.exportSolverView(), targetAmplitude, rightKd, curveStepCap, controlledIdx, errBound, derivatives, sys, conditionStopRcond);

    leftBranch = leftBranchTask.logs;
    rightBranch = rightBranchTask.logs;
    leftStatus = leftBranchTask.continuationStatus;
    rightStatus = rightBranchTask.continuationStatus;

    outputFile = fullfile(curveDataDir, ...
        sprintf('iso_amplitude_A%s_condition.mat', amplitude_tag(targetAmplitude)));
    save(outputFile, ...
        'targetAmplitude', ...
        'conditionStopRcond', ...
        'seed', ...
        'seedStatus', ...
        'leftBranch', ...
        'rightBranch', ...
        'leftStatus', ...
        'rightStatus');

    fprintf('Saved data: %s\n', outputFile);
end

%%
function task = solve_seed_task(baseSolverView, targetAmplitude, centerPeriod, controlledIdx, errBound, derivatives, sys, conditionStopRcond)
continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', conditionStopRcond);

itemsPer = struct([]);
itemsPer(1).prop = 'obsAmp';
itemsPer(1).idx = 1;
itemsPer(1).target = targetAmplitude;

itemsPer(2).prop = 'p_Psi';
itemsPer(2).idx = 1;
itemsPer(2).target = centerPeriod / (2 * pi);

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
function task = solve_branch_task(seedSolverView, targetAmplitude, targetKd, curveStepCap, controlledIdx, errBound, derivatives, sys, conditionStopRcond)
continuationOptions = struct( ...
    'initialLambdaStep', min(curveStepCap, 0.05), ...
    'maxLambdaStep', curveStepCap, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', conditionStopRcond);

itemsPer = struct([]);
itemsPer(1).prop = 'obsAmp';
itemsPer(1).idx = 1;
itemsPer(1).target = targetAmplitude;

itemsPer(2).prop = 'params';
itemsPer(2).idx = 1;
itemsPer(2).target = targetKd;

task = FMAM_ODE(sys, obs_spec(), seedSolverView, itemsPer, controlledIdx, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
task.psiUpdateMode = true;
task.refreshPsiModeReferences();
task.needLog = true;

task.fit()
task.step()
end

%%
function seedData = export_seed_data(task, y0)
finalParams = reshape(task.exportSolverView().params, 1, []);
orbitResult = circadian_forward_orbit(finalParams, y0, struct());
if ~orbitResult.success
    error('circadian_refactor:IsoAmplitudeSeedOrbitGenerationFailed', ...
        'Periodic-orbit extraction failed at the seed point.');
end

seedData = struct();
seedData.Parameters = finalParams;
seedData.period = orbitResult.features.period;
seedData.obsAmp = reshape(orbitResult.features.observable.amplitude, 1, []);
seedData.obsMax = reshape(orbitResult.features.observable.max, 1, []);
seedData.obsMin = reshape(orbitResult.features.observable.min, 1, []);
seedData.directConditionEstimate = task.logs(end).directConditionEstimate;
seedData.conditionNumber = 1 / seedData.directConditionEstimate;
end

%%
function tag = amplitude_tag(value)
tag = regexprep(sprintf('%.4f', value), '0+$', '');
tag = regexprep(tag, '\.$', '');
tag = strrep(tag, '.', 'p');
end

function obs = obs_spec()
obs = {@(variable) variable(:, 2) + variable(:, 3)};
end
