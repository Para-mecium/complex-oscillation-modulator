function result = circuit_forward_orbit(Parameters, y0, options)
if nargin < 3
    options = struct();
end

%% Default settings for the RLT circuit model
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));

sys = cell(1, 3);
run(fullfile(scriptDir, 'System.m'));

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
poOptions.solver_name = 'ode45';
poOptions.single_timespan = 500;
poOptions.max_windows = 3;
poOptions.event = 1;
poOptions.solver_tol = struct('RelTol', 1e-6, 'AbsTol', 1e-9);
poOptions.minCrossings = 6;
poOptions.transientFraction = 0.5;
poOptions.samplesPerCycle = 400;
poOptions.extractNumPoints = 500;

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

%% Forward integration and periodic-orbit extraction
rhs1 = sys{1};
rhs2 = sys{2};
rhs3 = sys{3};
odeFunc = @(t, y, parameter) [ ...
    rhs1(reshape(y, 1, []), parameter); ...
    rhs2(reshape(y, 1, []), parameter); ...
    rhs3(reshape(y, 1, []), parameter)];
poResult = extract_periodic_orbit(odeFunc, y0(:), reshape(Parameters, 1, []), poOptions);

result = struct();
result.Parameters = reshape(Parameters, 1, []);
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
