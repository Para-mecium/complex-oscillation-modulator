clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'data', 'common', 'initData.mat');
curveDataDir = fullfile(scriptDir, 'data', 'fig3c', 'curves');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(curveDataDir);

%% Base model and FMAM settings for Fig. 3c
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

%% Iso-period targets for Fig. 3c
targetPeriods = [50, 60, 70, 80, 90];
seedAmplitudes = [1, 2, 2.5, 2.5, 3];

leftI1 = 0.8;
rightI1 = 2.0;

curveStepCap = 2e-2;

%% Build system and derivatives
run(fullfile(scriptDir, 'System_base.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(baseInitialParams));

%% Build base FMAM state
baseState = state(obs, baseInitialParams, baseT, baseTSVar, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);
baseDerivedView = fmam_state_ops.derivedViewFromState(baseState);

%% Compute seed point and two end branches for each iso-period curve
for i = 1:numel(targetPeriods)
    outputFile = fullfile(curveDataDir, sprintf('fig3c_iso_period_T%03d.mat', round(targetPeriods(i))));

    % Generate seed
    continuationOptions = struct( ...
        'initialLambdaStep', 0.05,...
        'predictorMode', 'constant', ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', 1e-9);

    itemsPer = struct([]);
    itemsPer(1).prop = 'p_Psi';
    itemsPer(1).idx = 1;
    itemsPer(1).target = targetPeriods(i) / (2 * pi);

    itemsPer(2).prop = 'varAmp';
    itemsPer(2).idx = 2;
    itemsPer(2).target = seedAmplitudes(i);

    seedTask = FMAM_ODE(sys, obs, baseSolverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    seedTask.isPsiUpdated = true;
    seedTask.needLog = false;

    seedTask.fit()
    seedTask.step()
    seedTask.errBound = 1e-12;
    seedTask.fit()

    seed = struct();
    seed.params = reshape(seedTask.exportSolverView().params, 1, []);
    seed.solverView = seedTask.exportSolverView();
    seed.derivedView = seedTask.exportDerivedView();

    % left branch continuation
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

    % right branch continuation
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
