clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'data', 'common', 'initData.mat');
curveDataDir = fullfile(scriptDir, 'data', 'fig3d', 'curves');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(curveDataDir);

%% Base model and FMAM settings for Fig. 3d
initData = load(initDataFile);
baseInitialParams = initData.Parameters;
baseT = initData.TS{1};
baseTSVar = initData.TS{2};

obs = [];
M = 50;
PV.name = 'var';
PV.idx = 2;

controlledIdx = [1, 2];
errBound = 1e-6;

%% Iso-amplitude targets for Fig. 3d
targetAmplitudes = 1.2:0.3:3.6;
seedPeriods = [50:3:71, 80];

leftI1 = 0.8;
rightI1 = 2.0;

curveStepCap = 5e-3;

maxParamJump = [0.2, 0.25];
maxPeriodJump = 15;

%% Build system and derivatives
run(fullfile(scriptDir, 'System_base.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(baseInitialParams));

%% Build base FMAM state
baseState = state(obs, baseInitialParams, baseT, baseTSVar, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);

%% Compute seed point and two end branches for each iso-amplitude curve
for i = 1:numel(targetAmplitudes)
    amplitudeTag = strrep(sprintf('%.1f', targetAmplitudes(i)), '.', 'p');
    outputFile = fullfile(curveDataDir, sprintf('fig3d_iso_amplitude_A%s.mat', amplitudeTag));

    continuationOptions = struct( ...
        'initialLambdaStep', 0.05, ...
        'predictorMode', 'constant', ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', 1e-9);

    itemsPer = struct([]);
    itemsPer(1).prop = 'varAmp';
    itemsPer(1).idx = 2;
    itemsPer(1).target = targetAmplitudes(i);

    itemsPer(2).prop = 'p_Psi';
    itemsPer(2).idx = 1;
    itemsPer(2).target = seedPeriods(i) / (2 * pi);

    seedTask = FMAM_ODE(sys, obs, baseSolverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    seedTask.isPsiUpdated = true;
    seedTask.needLog = false;

    seedTask.fit()
    seedTask.step()
    seedTask.errBound = 1e-12;
    seedTask.fit()

    % left branch
    seed = struct();
    seed.params = reshape(seedTask.exportSolverView().params, 1, []);
    seed.solverView = seedTask.exportSolverView();
    seed.derivedView = seedTask.exportDerivedView();

    continuationOptions = struct( ...
        'predictorMode', 'constant', ...
        'initialLambdaStep', min(curveStepCap, 0.05), ...
        'maxLambdaStep', curveStepCap, ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', 1e-9);

    itemsPer(2).prop = 'params';
    itemsPer(2).idx = 1;
    itemsPer(2).target = leftI1;

    leftBranchTask = FMAM_ODE(sys, obs, seed.solverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    leftBranchTask.isPsiUpdated = true;
    leftBranchTask.needLog = true;

    leftBranchTask.fit()
    leftBranchTask.step()

    leftBranch = leftBranchTask.logs;

    % right branch

    continuationOptions = struct( ...
        'predictorMode', 'constant', ...
        'initialLambdaStep', min(curveStepCap, 0.05), ...
        'maxLambdaStep', curveStepCap, ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', 1e-9);

    itemsPer(2).prop = 'params';
    itemsPer(2).idx = 1;
    itemsPer(2).target = rightI1;

    rightBranchTask = FMAM_ODE(sys, obs, seed.solverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    rightBranchTask.isPsiUpdated = true;
    rightBranchTask.needLog = true;

    rightBranchTask.fit()
    rightBranchTask.step()

    rightBranch = rightBranchTask.logs;

    save(outputFile, 'seed', 'leftBranch', 'rightBranch');
    fprintf('Saved data: %s\n', outputFile);
end