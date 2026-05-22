clear
clc

%% Local configuration
% evalBudget is the total number of calls to circuit_forward_orbit.
% If refinement is enabled, the search budget is evalBudget - refinementBudget.
evalBudget = 1000;
refinementBudget = 300;
refinementTopK = 2;
randomSeeds = 1:100;
numWorkers = 8;
methods = { ...
    struct('name', 'uniform', 'refinement', 'powell'), ...
    struct('name', 'lhs', 'refinement', 'powell'), ...
    struct('name', 'de', 'refinement', 'none'), ...
    struct('name', 'sobol', 'refinement', 'powell'), ...
    struct('name', 'pso', 'refinement', 'none')};

lossName = 'relative_l2_orbit'; % property_difference, relative_l2_orbit
penaltyLoss = 1e30;
saveBestOrbit = true;
odeSolverTolerance = struct('RelTol', 1e-6, 'AbsTol', 1e-6);
progressEnabled = true;
progressInterval = 50;

%% Shared setup
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');
addpath(parameterInferenceDir, '-begin');
addpath(scriptDir, '-begin');

userConfig = struct();
userConfig.evalBudget = evalBudget;
userConfig.randomSeeds = randomSeeds;
userConfig.refinement.budget = refinementBudget;
userConfig.refinement.topK = refinementTopK;
userConfig.lossName = lossName;
userConfig.penaltyLoss = penaltyLoss;
userConfig.saveBestOrbit = saveBestOrbit;
userConfig.forwardOptions.poOptions.solver_tol = odeSolverTolerance;
userConfig.progress.enabled = progressEnabled;
userConfig.progress.interval = progressInterval;
userConfig.targetDataFile = fullfile(parameterInferenceDir, 'initData_circuit.mat');
userConfig.baseDataFile = fullfile(parameterInferenceDir, 'initData_ODE.mat');
userConfig.resultDir = fullfile(scriptDir, 'results');

config = build_config(userConfig);

if ~isfolder(config.resultDir)
    mkdir(config.resultDir);
end

summaryTemplate = struct( ...
    'methodName', '', ...
    'refinementMethod', '', ...
    'seed', NaN, ...
    'resultFile', '', ...
    'bestLoss', NaN, ...
    'bestEvalIndex', NaN, ...
    'numEvaluations', NaN, ...
    'numSuccessfulForward', NaN, ...
    'refinementEnabled', false, ...
    'refinedBestLoss', NaN, ...
    'refinementEvaluations', 0, ...
    'totalEvaluations', NaN, ...
    'runtime', NaN, ...
    'totalRuntime', NaN, ...
    'status', '', ...
    'message', '');
tasks = build_tasks(methods, randomSeeds);
summaryRows = repmat(summaryTemplate, numel(tasks), 1);

fprintf(['Baseline comparison starts: %d tasks (%d methods x %d seeds), ' ...
    'budget=%d, workers=%d.\n'], ...
    numel(tasks), numel(methods), numel(randomSeeds), config.evalBudget, numWorkers);

if numWorkers <= 1
    for taskIdx = 1:numel(tasks)
        summaryRows(taskIdx) = run_one_task(tasks(taskIdx), config, summaryTemplate, false);
    end
else
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= numWorkers
        if ~isempty(pool)
            delete(pool);
        end
        parpool(numWorkers);
    end

    parfor taskIdx = 1:numel(tasks)
        summaryRows(taskIdx) = run_one_task(tasks(taskIdx), config, summaryTemplate, true);
    end
end

summaryFile = fullfile(config.resultDir, ...
    sprintf('baseline_comparison_summary_%s.mat', result_name_token(config.lossName)));
configSummary = config_summary(config);
save(summaryFile, 'summaryRows', 'configSummary');
fprintf('Saved baseline comparison summary: %s\n', summaryFile);

function tasks = build_tasks(methods, randomSeeds)
numTasks = numel(methods) * numel(randomSeeds);
tasks = repmat(struct('methodSpec', struct(), 'seed', NaN), numTasks, 1);
taskIdx = 0;
for seedIdx = 1:numel(randomSeeds)
    for methodIdx = 1:numel(methods)
        taskIdx = taskIdx + 1;
        tasks(taskIdx).methodSpec = methods{methodIdx};
        tasks(taskIdx).seed = randomSeeds(seedIdx);
    end
end
end

function row = run_one_task(task, config, summaryTemplate, disableProgress)
addpath(config.repoDir, '-begin');
addpath(config.circuitDir, '-begin');
addpath(config.parameterInferenceDir, '-begin');
addpath(config.scriptDir, '-begin');

methodSpec = task.methodSpec;
seed = task.seed;
methodName = methodSpec.name;
refinementMethod = get_refinement_method(methodSpec);
useRefinement = ~strcmp(refinementMethod, 'none');

if useRefinement
    if config.refinement.budget >= config.evalBudget
        error('run_baseline_comparison:InvalidBudget', ...
            ['refinementBudget must be smaller than evalBudget when ' ...
            'refinement is enabled.']);
    end
    methodConfig = config;
    methodConfig.evalBudget = config.evalBudget - config.refinement.budget;
else
    methodConfig = config;
end

if disableProgress
    methodConfig.progress.enabled = false;
end

fprintf(['[%s] seed=%d starts: searchBudget=%d, ' ...
    'refinementBudget=%d, totalBudget=%d.\n'], ...
    methodName, seed, methodConfig.evalBudget, ...
    useRefinement * config.refinement.budget, config.evalBudget);

row = summaryTemplate;
row.methodName = methodName;
row.refinementMethod = refinementMethod;
row.seed = seed;

try
    result = run_search_method(methodName, methodConfig, seed);

    result.totalEvalBudget = config.evalBudget;
    result.searchEvalBudget = methodConfig.evalBudget;
    result.refinementMethod = refinementMethod;
    result.searchRuntime = result.runtime;
    if useRefinement
        result.refinementBudget = config.refinement.budget;
    else
        result.refinementBudget = 0;
    end

    refinementConfig = config;
    if disableProgress
        refinementConfig.progress.enabled = false;
    end

    result.refinement = run_refinement(refinementMethod, refinementConfig, seed, result);
    result.refinementRuntime = result.refinement.runtime;
    result.runtime = result.searchRuntime + result.refinementRuntime;
    result.totalBestSoFarLosses = build_total_best_so_far(result);
    result.totalTimeSoFar = build_total_time_so_far(result);

    methodResultDir = result_method_loss_dir(config.resultDir, methodName, config.lossName);
    if ~isfolder(methodResultDir)
        mkdir(methodResultDir);
    end

    resultFile = fullfile(methodResultDir, ...
        sprintf('baseline_%s_refine_%s_seed_%d.mat', ...
        methodName, refinementMethod, seed));
    configSummary = config_summary(config);
    save(resultFile, 'result', 'configSummary');

    row.resultFile = resultFile;
    row.bestLoss = result.bestLoss;
    row.bestEvalIndex = result.bestEvalIndex;
    row.numEvaluations = result.numEvaluations;
    row.numSuccessfulForward = nnz(result.successFlags);
    row.refinementEnabled = result.refinement.enabled;
    row.refinedBestLoss = result.refinement.bestLoss;
    row.refinementEvaluations = result.refinement.numEvaluations;
    row.totalEvaluations = result.numEvaluations + result.refinement.numEvaluations;
    row.runtime = result.runtime;
    row.totalRuntime = result.runtime;
    row.status = 'success';
    row.message = '';

    fprintf(['[%s] seed=%d done: bestLoss=%.6e, refinedBestLoss=%.6e, ' ...
        'eval=%d, totalEval=%d, saved=%s\n'], ...
        methodName, seed, result.bestLoss, row.refinedBestLoss, ...
        result.bestEvalIndex, row.totalEvaluations, resultFile);
catch ME
    row.status = 'failed';
    row.message = ME.message;
    fprintf('[%s] seed=%d failed: %s\n', methodName, seed, ME.message);
end
end

function result = run_search_method(methodName, config, seed)
switch methodName
    case 'uniform'
        result = baseline_uniform(config, seed);
    case 'lhs'
        result = baseline_lhs(config, seed);
    case 'de'
        result = baseline_de(config, seed);
    case 'sobol'
        result = baseline_sobol(config, seed);
    case 'pso'
        result = baseline_pso(config, seed);
    otherwise
        error('run_baseline_comparison:UnknownMethod', ...
            'Unknown baseline method: %s.', methodName);
end
end

function refinementMethod = get_refinement_method(methodSpec)
if ~isfield(methodSpec, 'refinement') || isempty(methodSpec.refinement)
    refinementMethod = 'none';
    return
end

refinementMethod = lower(string(methodSpec.refinement));
refinementMethod = char(refinementMethod);
if ~ismember(refinementMethod, {'none', 'powell', 'fmincon'})
    error('run_baseline_comparison:UnknownRefinement', ...
        'Unknown refinement method: %s.', refinementMethod);
end
end

function refinement = run_refinement(refinementMethod, config, seed, result)
switch refinementMethod
    case 'none'
        refinement = disabled_refinement(config);
    case 'powell'
        refinement = run_powell(config, seed, result);
    case 'fmincon'
        refinement = run_fmincon(config, seed, result);
    otherwise
        error('run_baseline_comparison:UnknownRefinement', ...
            'Unknown refinement method: %s.', refinementMethod);
end
end

function totalBestSoFarLosses = build_total_best_so_far(result)
if result.refinement.numEvaluations > 0
    totalBestSoFarLosses = [ ...
        result.bestSoFarLosses(:); ...
        result.refinement.globalBestSoFarLosses(:)];
else
    totalBestSoFarLosses = result.bestSoFarLosses(:);
end
end

function totalTimeSoFar = build_total_time_so_far(result)
searchTime = result.timeSoFar(:);
if result.refinement.numEvaluations > 0
    totalTimeSoFar = [ ...
        searchTime; ...
        result.searchRuntime + result.refinement.timeSoFar(:)];
else
    totalTimeSoFar = searchTime;
end
end
