function result = HH_forward_orbit(p, I, y0, N, options)
if nargin < 5
    options = struct();
end

%% Default settings for the HH model
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));

obs = {};
observableSpec = [];
observableFile = fullfile(scriptDir, 'Observable.m');
if isfile(observableFile)
    run(observableFile);
    if exist('obs', 'var') && ~isempty(obs)
        observableSpec = @(t, y) obs{1}(y);
    end
end

poOptions = struct();
poOptions.orbit_solver = 'direct';
poOptions.solver_name = 'ode15s';
poOptions.single_timespan = 1000;
poOptions.max_windows = 10;
poOptions.event = 1;
poOptions.solver_tol = struct('RelTol', 1e-6, 'AbsTol', 1e-8);
poOptions.minCrossings = 6;
poOptions.transientFraction = 0.5;
poOptions.samplesPerCycle = 500;
poOptions.extractNumPoints = 600;

featureOptions = struct();
featureOptions.extremaRefinement = true;

%% Override defaults only where this interface needs it
if isfield(options, 'poOptions') && ~isempty(options.poOptions)
    optionNames = fieldnames(options.poOptions);
    for i = 1:numel(optionNames)
        name = optionNames{i};
        poOptions.(name) = options.poOptions.(name);
    end
end

if isfield(options, 'featureOptions') && ~isempty(options.featureOptions)
    optionNames = fieldnames(options.featureOptions);
    for i = 1:numel(optionNames)
        name = optionNames{i};
        featureOptions.(name) = options.featureOptions.(name);
    end
end

%% Forward integration and periodic-orbit extraction
parameter = struct( ...
    'p', p, ...
    'I', I(:), ...
    'N', N);

odeFunc = @(t, y, parameterStruct) build_model( ...
    t, y, parameterStruct.I, parameterStruct.p, parameterStruct.N);
poResult = extract_periodic_orbit(odeFunc, y0(:), parameter, poOptions);

result = struct();
result.p = p;
result.I = I(:);
result.N = N;
result.initialState = y0(:);
result.poResult = poResult;

if ~poResult.has_orbit
    result.success = false;
    result.msg = poResult.message;
    return
end

orbit = struct();
orbit.t = poResult.orbit_t(:);
orbit.y = poResult.orbit_y;
orbit.period = poResult.period;

features = evaluate_orbit_features(orbit, observableSpec, [], featureOptions);

result.success = true;
result.msg = "";
result.orbit = orbit;
result.features = features;
end
