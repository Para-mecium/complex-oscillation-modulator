function config = build_config(userConfig)
if nargin < 1
    userConfig = struct();
end

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');
addpath(scriptDir, '-begin');

config = struct();
config.scriptDir = scriptDir;
config.parameterInferenceDir = parameterInferenceDir;
config.circuitDir = circuitDir;
config.repoDir = repoDir;

config.evalBudget = 1000;
config.randomSeeds = 0:9;
config.targetDataFile = fullfile(parameterInferenceDir, 'initData_ODE.mat');
config.baseDataFile = fullfile(parameterInferenceDir, 'initData_ODE.mat');
config.resultDir = fullfile(scriptDir, 'results');

config.activeIndex = [1 4 5 6];
config.activeLowerBound = [500 0.5 0.5 0.5];
config.activeUpperBound = [1500 2 2 2];
config.activeNames = {'R_C', 'inv_C_1', 'inv_C_2', 'inv_C_3'};

config.lossName = 'relative_l2_orbit';
config.lossOptions = struct();
config.lossOptions.name = config.lossName;
config.lossOptions.compareNumPoints = 500;
config.lossOptions.periodWeight = 0;
config.lossOptions.phaseAlignment = true;

config.penaltyLoss = 1e30;
config.forwardOptions = struct();
config.forwardOptions.poOptions = struct();
config.forwardOptions.poOptions.solver_tol = struct('RelTol', 1e-6, 'AbsTol', 1e-9);
config.saveBestOrbit = true;

config.refinement = struct();
config.refinement.budget = 200;
config.refinement.topK = 3;
config.refinement.maxSweeps = 20;
config.refinement.lineSearchMaxEval = 12;
config.refinement.tolerance = 1e-4;
config.refinement.lossTolerance = 1e-8;
config.refinement.minLineSearchWidth = 1e-4;
config.refinement.fminconAlgorithm = 'interior-point';
config.refinement.fminconStepTolerance = 1e-6;
config.refinement.fminconOptimalityTolerance = 1e-6;

config.progress = struct();
config.progress.enabled = true;
config.progress.interval = [];

config.de = struct();
config.de.populationSize = [];
config.de.F = 0.8;
config.de.CR = 0.9;

config.pso = struct();
config.pso.swarmSize = [];
config.pso.inertia = 0.72;
config.pso.cognitiveWeight = 1.49;
config.pso.socialWeight = 1.49;
config.pso.initialVelocityScale = 0.1;

config = merge_user_config(config, userConfig);
config.lossOptions.name = config.lossName;
if isempty(config.progress.interval)
    config.progress.interval = max(1, floor(config.evalBudget / 20));
end

baseData = load(config.baseDataFile);
config.baseParameters = reshape(baseData.Parameters, 1, []);
config.y0 = baseData.TS{2}(1, :).';

targetData = load(config.targetDataFile);
if isfield(targetData, 't') && isfield(targetData, 'y')
    config.targetOrbit = struct( ...
        't', targetData.t(:), ...
        'y', targetData.y, ...
        'period', read_target_period(targetData));
    config.targetFeatures = read_target_features(targetData);
    config.targetParameters = [];
elseif isfield(targetData, 'TS') && isfield(targetData, 'Parameters')
    targetParameters = reshape(targetData.Parameters, 1, []);
    targetInitialState = targetData.TS{2}(1, :).';
    targetForward = circuit_forward_orbit( ...
        targetParameters, targetInitialState, config.forwardOptions);
    if ~targetForward.success
        error('build_config:TargetOrbitFailed', ...
            'Target periodic-orbit extraction failed: %s', targetForward.msg);
    end

    config.targetOrbit = targetForward.orbit;
    config.targetFeatures = targetForward.features;
    config.targetParameters = targetParameters;
else
    error('build_config:UnsupportedTargetData', ...
        'Unsupported target data format: %s', config.targetDataFile);
end
end

function config = merge_user_config(config, userConfig)
names = fieldnames(userConfig);
for i = 1:numel(names)
    name = names{i};
    if isstruct(userConfig.(name)) && isfield(config, name) && isstruct(config.(name))
        nestedNames = fieldnames(userConfig.(name));
        for j = 1:numel(nestedNames)
            nestedName = nestedNames{j};
            config.(name).(nestedName) = userConfig.(name).(nestedName);
        end
    else
        config.(name) = userConfig.(name);
    end
end
end

function period = read_target_period(targetData)
if isfield(targetData, 'period') && ~isempty(targetData.period)
    period = targetData.period;
else
    t = targetData.t(:);
    period = t(end) - t(1);
end
end

function features = read_target_features(targetData)
features = struct();
if isfield(targetData, 'period')
    features.period = targetData.period;
end
if isfield(targetData, 'varAmp') || isfield(targetData, 'varMax') || isfield(targetData, 'varMin')
    features.state = struct();
    if isfield(targetData, 'varAmp')
        features.state.amplitude = targetData.varAmp;
    end
    if isfield(targetData, 'varMax')
        features.state.max = targetData.varMax;
    end
    if isfield(targetData, 'varMin')
        features.state.min = targetData.varMin;
    end
end
end
