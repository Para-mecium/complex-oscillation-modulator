function result = flexmod_forward_orbit(activeParams, initialState, options)

%% Default settings for the flexible modulator model
scriptDir = fileparts(mfilename('fullpath'));
sys = cell(1, 2);
if strcmp(options.systemName, 'base')
    run(fullfile(scriptDir, 'System_base.m'));
elseif strcmp(options.systemName, 'temp')
    run(fullfile(scriptDir, 'System_temp.m'));
end

poOptions = struct();
poOptions.orbit_solver = 'direct';
poOptions.solver_name = 'ode45';
poOptions.single_timespan = 1600;
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
odeFunc = @(t, y, parameter) [ ...
    rhs1(reshape(y, 1, []), parameter); ...
    rhs2(reshape(y, 1, []), parameter)];
poResult = extract_periodic_orbit(odeFunc, initialState, activeParams, poOptions);

result = struct();
result.activeParams = reshape(activeParams, 1, []);
result.initialState = initialState(:);
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

features = evaluate_orbit_features(orbit, [], [], featureOptions);

result.success = true;
result.orbit = orbit;
result.features = features;
end
