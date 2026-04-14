function result = Modulation_bi_target(varargin)
%MODULATION_BI_TARGET Build self-contained bi-target modulation caches.
%   MODULATION_BI_TARGET() generates reusable periodic-orbit caches for
%   the 1.0x, 1.5x, and 2.0x period targets while preserving the V1
%   amplitude from the learned baseline in learnedData_ODE.mat.
%
%   MODULATION_BI_TARGET('generatePathCache', true) also generates the
%   dense continuation-path cache used by draw_params.m.

opts = struct();
opts.generatePathCache = false;
opts.lambdaStepCap = 5e-3;
for i = 1:2:numel(varargin)
    name = varargin{i};
    value = varargin{i + 1};
    switch lower(name)
        case 'generatepathcache'
            opts.generatePathCache = logical(value);
        case 'lambdastepcap'
            opts.lambdaStepCap = double(value);
        otherwise
            error('Modulation_bi_target:UnknownOption', ...
                'Unknown option %s.', char(string(name)));
    end
end

scriptDir = fileparts(mfilename('fullpath'));
circuitDir = fileparts(scriptDir);
fmamDir = fileparts(circuitDir);

startDir = pwd;
cleanupObj = onCleanup(@() cd(startDir));

rawData = load(fullfile(circuitDir, 'learnedData_ODE.mat'));

targetMultipliers = [1, 1.5, 2];
targetFiles = { ...
    fullfile(scriptDir, 'period_target_1p0x.mat'), ...
    fullfile(scriptDir, 'period_target_1p5x.mat'), ...
    fullfile(scriptDir, 'period_target_2p0x.mat')};

pathFile = fullfile(scriptDir, 'params_modulation_path.mat');
needPath = opts.generatePathCache;
M = 75;
PV = struct('name', 'var', 'idx', 1);
params = reshape(double(rawData.Parameters), 1, []);
t = rawData.TS{1};
TS_var = rawData.TS{2};
obs = [];

cd(circuitDir);
System

cd(fmamDir);
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

baseState = state(obs, params, t, TS_var, M, PV);
baseState.updatePeriod();
baseState.updateVar2();

ctx = struct();
ctx.sys = sys;
ctx.obs = obs;
ctx.derivatives = derivatives;
ctx.params = params;
ctx.t = t;
ctx.TS_var = TS_var;
ctx.M = M;
ctx.PV = PV;
ctx.basePeriod = baseState.period;
ctx.baseVarAmp1 = double(baseState.varAmp(1));
ctx.poExtractDir = fullfile(fmamDir, 'PO_extract');

Parameters = reshape(double(rawData.Parameters), 1, []);
TS = rawData.TS;
period = double(rawData.period);
varAmp = reshape(double(rawData.varAmp), 1, []);
varMax = reshape(double(rawData.varMax), 1, []);
varMin = reshape(double(rawData.varMin), 1, []);
periodMultiplier = 1;
targetPeriod = period;
save(targetFiles{1}, ...
    'Parameters', 'TS', 'period', 'varAmp', 'varMax', 'varMin', ...
    'periodMultiplier', 'targetPeriod');

for i = 2:numel(targetMultipliers)
    periodMultiplier = targetMultipliers(i);
    stateObj = create_state(ctx);
    [items_per, items_controlled, ~, periodTarget] = ...
        build_bi_target_setup(stateObj, ctx.basePeriod, ctx.baseVarAmp1, periodMultiplier);

    errBound = 1e-6;
    continuationOptions = struct('initialLambdaStep', 0.01);

    task = FMAM_ODE(ctx.sys, ctx.obs, stateObj, items_per, items_controlled, ...
        [], errBound, 'derivatives', ctx.derivatives, ...
        'continuationOptions', continuationOptions);
    task.isPsiUpdated = true;

    task.fit();
    task.step();
    task.errBound = 1e-12;
    task.fit();

    stateView = task.exportSolverView();
    stateDerived = task.exportDerivedView();

    Parameters = reshape(double(stateView.params), 1, []);
    addpath(ctx.poExtractDir);

    odeFunc = @(t, y, parameter) [ ...
        ctx.sys{1}(y(:).', parameter); ...
        ctx.sys{2}(y(:).', parameter); ...
        ctx.sys{3}(y(:).', parameter)];
    y0 = stateDerived.TS_var(1, :).';
    searchWindow = max(10 * stateDerived.period, stateDerived.t(end));
    orbitOptions = struct( ...
        'solver_name', 'ode45', ...
        'tspan', [0, searchWindow], ...
        'event', 1, ...
        'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-6), ...
        'minCrossings', 3, ...
        'transientFraction', 0);

    poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
    if ~poResult.has_orbit
        error('Modulation_bi_target:PeriodicOrbitGenerationFailed', ...
            'Periodic-orbit extraction failed for %.1fx (%s).', ...
            periodMultiplier, poResult.message);
    end

    poState = state(ctx.obs, Parameters, poResult.orbit_t, poResult.orbit_y, ctx.M, ctx.PV);
    poState.updatePeriod();
    poState.updateVar2();

    TS = {poState.t, poState.TS_var};
    period = poState.period;
    varAmp = reshape(poState.varAmp, 1, []);
    varMax = reshape(poState.varMax, 1, []);
    varMin = reshape(poState.varMin, 1, []);
    save(targetFiles{i}, ...
        'Parameters', 'TS', 'period', 'varAmp', 'varMax', 'varMin', ...
        'periodMultiplier', 'periodTarget');
end

if needPath
    stateObj = create_state(ctx);
    [items_per, items_controlled, ~] = build_bi_target_setup( ...
        stateObj, ctx.basePeriod, ctx.baseVarAmp1, 2);

    errBound = 1e-6;
    continuationOptions = struct( ...
        'initialLambdaStep', opts.lambdaStepCap, ...
        'maxLambdaStep', opts.lambdaStepCap);

    task = FMAM_ODE(ctx.sys, ctx.obs, stateObj, items_per, items_controlled, ...
        [], errBound, 'derivatives', ctx.derivatives, ...
        'continuationOptions', continuationOptions);
    task.isPsiUpdated = true;
    task.needLog = true;

    task.fit();
    params_start = reshape(double(stateObj.params), 1, []);

    task.step();
    task.errBound = 1e-12;
    task.fit();

    stateView = task.exportSolverView();
    params_end = reshape(double(stateView.params), 1, []);

    solution = task.logs;
    if isempty(solution)
        error('Modulation_bi_target:MissingContinuationLogs', ...
            'No continuation logs were captured for the path cache.');
    end

    curve_params = zeros(numel(solution), numel(solution(1).params));
    for i = 1:numel(solution)
        curve_params(i, :) = reshape(double(solution(i).params), 1, []);
    end

    pathData = struct();
    pathData.curve_params = curve_params;
    pathData.params_start = params_start;
    pathData.params_end = params_end;
    save(pathFile, '-struct', 'pathData');
end

result = struct();
result.targetFiles = targetFiles;
result.pathFile = pathFile;
end

function stateObj = create_state(ctx)
stateObj = state(ctx.obs, ctx.params, ctx.t, ctx.TS_var, ctx.M, ctx.PV);
stateObj.updatePeriod();
stateObj.updateVar2();
end

function [items_per, items_controlled, max_stepsize, periodTarget] = ...
        build_bi_target_setup(stateObj, basePeriod, baseVarAmp1, periodMultiplier)
periodTarget = basePeriod * periodMultiplier;

items_per = struct([]);
items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = periodTarget / (2 * pi);

items_per(2).prop = 'varAmp';
items_per(2).idx = 1;
items_per(2).target = baseVarAmp1;

items_controlled = [1, 4];

stepsNum = 100;
max_stepsize_const = [1e-3, 1e-2];
targetDelta = [ ...
    (periodTarget - stateObj.period) / (2 * pi), ...
    baseVarAmp1 - stateObj.varAmp(1)];
max_stepsize = max(abs(targetDelta) / stepsNum, max_stepsize_const);
end
