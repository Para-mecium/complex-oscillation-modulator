clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
flexmodDir = fileparts(scriptDir);
repoDir = fileparts(flexmodDir);
initDataFile = fullfile(flexmodDir, 'data', 'common', 'initData.mat');
curveDataDir = fullfile(scriptDir, 'data', 'iso_period', 'curves');

addpath(repoDir);
addpath(flexmodDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(curveDataDir);

%% Base model and FMAM settings
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

conditionStopRcond = 1e-12;

%% Iso-period targets
targetPeriods = [50, 60, 70, 80, 90];
seedAmplitudes = [1, 2, 2.5, 2.5, 3];

leftI1 = 0.8;
rightI1 = 2.0;

curveStepCap = 2e-2;

%% Build system and derivatives
run(fullfile(flexmodDir, 'System_base.m'));
derivatives = build_symbolic_derivatives(sys, obs, numel(baseInitialParams));

%% Build base FMAM state
baseState = state(obs, baseInitialParams, baseT, baseTSVar, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);

%% Compute seed point and two end branches for each iso-period curve
for i = 1:numel(targetPeriods)
    outputFile = fullfile(curveDataDir, sprintf('iso_period_T%03d_condition.mat', round(targetPeriods(i))));

    continuationOptions = struct( ...
        'initialLambdaStep', 0.05,...
        'predictorMode', 'constant', ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', conditionStopRcond);

    itemsPer = struct([]);
    itemsPer(1).prop = 'p_Psi';
    itemsPer(1).idx = 1;
    itemsPer(1).target = targetPeriods(i) / (2 * pi);

    itemsPer(2).prop = 'varAmp';
    itemsPer(2).idx = 2;
    itemsPer(2).target = seedAmplitudes(i);

    seedTask = FMAM_ODE(sys, obs, baseSolverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    seedTask.psiUpdateMode = true;
    seedTask.refreshPsiModeReferences();
    seedTask.needLog = true;

    seedTask.fit()
    seedTask.step()
    seedTask.errBound = 1e-12;
    seedFit = seedTask.fit();

    seed = struct();
    seed.params = reshape(seedTask.exportSolverView().params, 1, []);
    seed.solverView = seedTask.exportSolverView();
    seed.derivedView = seedTask.exportDerivedView();
    seed.directConditionEstimate = seedTask.logs(end).directConditionEstimate;
    seed.conditionNumber = 1 / seed.directConditionEstimate;

    continuationOptions = struct( ...
        'predictorMode', 'constant', ...
        'initialLambdaStep', min(curveStepCap, 0.05), ...
        'maxLambdaStep', curveStepCap, ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', conditionStopRcond);

    itemsPer(2).prop = 'params';
    itemsPer(2).idx = 1;
    itemsPer(2).target = leftI1;

    leftBranchTask = FMAM_ODE(sys, obs, seed.solverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    leftBranchTask.psiUpdateMode = true;
    leftBranchTask.refreshPsiModeReferences();
    leftBranchTask.needLog = true;

    leftBranchTask.fit()
    leftBranchTask.step()

    leftBranch = leftBranchTask.logs;
    leftStatus = leftBranchTask.continuationStatus;

    continuationOptions = struct( ...
        'predictorMode', 'constant', ...
        'initialLambdaStep', min(curveStepCap, 0.05), ...
        'maxLambdaStep', curveStepCap, ...
        'conditionStopEnabled', true, ...
        'conditionStopRcond', conditionStopRcond);

    itemsPer(2).prop = 'params';
    itemsPer(2).idx = 1;
    itemsPer(2).target = rightI1;

    rightBranchTask = FMAM_ODE(sys, obs, seed.solverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    rightBranchTask.psiUpdateMode = true;
    rightBranchTask.refreshPsiModeReferences();
    rightBranchTask.needLog = true;

    rightBranchTask.fit()
    rightBranchTask.step()

    rightBranch = rightBranchTask.logs;
    rightStatus = rightBranchTask.continuationStatus;

    save(outputFile, ...
        'targetPeriods', ...
        'seedAmplitudes', ...
        'conditionStopRcond', ...
        'seed', ...
        'leftBranch', ...
        'rightBranch', ...
        'leftStatus', ...
        'rightStatus');
    fprintf('Saved data: %s\n', outputFile);
end
