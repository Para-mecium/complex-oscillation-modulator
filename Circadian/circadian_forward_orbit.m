function result = circadian_forward_orbit(Parameters, y0, options)
if nargin < 3
    options = struct();
end

%% Default settings for the circadian model
scriptDir = fileparts(mfilename('fullpath'));
sys = cell(1, 3);
run(fullfile(scriptDir, 'System.m'));

poOptions = struct();
poOptions.orbit_solver = 'direct';
poOptions.solver_name = 'ode45';
poOptions.single_timespan = 400;
poOptions.max_windows = 3;
poOptions.event = 1;
poOptions.solver_tol = struct('RelTol', 1e-6, 'AbsTol', 1e-9);
poOptions.minCrossings = 6;
poOptions.transientFraction = 0.5;
poOptions.samplesPerCycle = 400;
poOptions.extractNumPoints = 500;
poOptions.observableFcn = @observable_ptot;

featureOptions = struct();
featureOptions.extremaRefinement = true;

%% Override periodic-orbit extraction settings when explicitly provided
if isfield(options, 'poOptions') && ~isempty(options.poOptions)
    optionNames = fieldnames(options.poOptions);
    for i = 1:numel(optionNames)
        name = optionNames{i};
        poOptions.(name) = options.poOptions.(name);
    end
end

function value = observable_ptot(~, y, ~, ~)
value = y(2) + y(3);
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

if ~poResult.has_orbit
    result.success = false;
    result.message = poResult.message;
    return
end

orbit = struct();
orbit.t = poResult.orbit_t(:);
orbit.y = poResult.orbit_y;
orbit.period = poResult.period;

observableSpec = @(t, y) y(:, 2) + y(:, 3); 
featureSpec = ["period", "observable_max", "observable_min", "observable_amplitude"];
features = evaluate_orbit_features(orbit, observableSpec, featureSpec, featureOptions);

result.success = true;
result.orbit = orbit;
result.features = features;
end
