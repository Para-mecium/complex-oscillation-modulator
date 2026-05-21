function result = repressilator_forward_orbit(Parameters, y0, options)
if nargin < 3
    options = struct();
end

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));

sys = cell(1, 6);
run(fullfile(scriptDir, 'System.m'));

poOptions = struct();
poOptions.orbit_solver = 'direct';
poOptions.solver_name = 'ode45';
poOptions.tspan = [0, 2000];
poOptions.event = 4;
poOptions.solver_tol = struct('RelTol', 1e-8, 'AbsTol', 1e-10);
poOptions.minCrossings = 6;
poOptions.transientFraction = 0.5;
poOptions.samplesPerCycle = 400;
poOptions.extractNumPoints = 500;

if isfield(options, 'poOptions') && ~isempty(options.poOptions)
    names = fieldnames(options.poOptions);
    for i = 1:numel(names)
        poOptions.(names{i}) = options.poOptions.(names{i});
    end
end

odeFunc = @(t, y, parameter) [ ...
    sys{1}(reshape(y, 1, []), parameter); ...
    sys{2}(reshape(y, 1, []), parameter); ...
    sys{3}(reshape(y, 1, []), parameter); ...
    sys{4}(reshape(y, 1, []), parameter); ...
    sys{5}(reshape(y, 1, []), parameter); ...
    sys{6}(reshape(y, 1, []), parameter)];
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
features = evaluate_orbit_features(orbit, [], [], struct());

result.success = true;
result.msg = "";
result.orbit = orbit;
result.features = features;
end
