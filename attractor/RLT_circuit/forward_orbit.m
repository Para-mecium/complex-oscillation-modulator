function result = forward_orbit(Parameters, y0, options)
%FORWARD_ORBIT Extract a periodic orbit for the RLT circuit model.

if nargin < 1 || isempty(Parameters)
    Parameters = [500, 6.6, 0.06, 1, 1, 1, 0.00297];
end
if nargin < 2 || isempty(y0)
    y0 = [1.383080070970441; 0.07615547480017665; 0.04041540879836707];
end
if nargin < 3 || isempty(options)
    options = struct();
end

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(fileparts(scriptDir));
addpath(scriptDir, '-begin');
addpath(fullfile(repoDir, 'PO_extract'), '-begin');

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
poOptions = merge_struct(poOptions, get_option_struct(options, 'poOptions'));

featureOptions = struct();
featureOptions.extremaRefinement = true;
featureOptions = merge_struct(featureOptions, get_option_struct(options, 'featureOptions'));

Parameters = reshape(Parameters, 1, []);
y0 = y0(:);
odeFunc = @(t, y, parameter) build_model(t, y, parameter);
poResult = extract_periodic_orbit(odeFunc, y0, Parameters, poOptions);

result = struct();
result.Parameters = Parameters;
result.initialState = y0;
result.poResult = poResult;

if ~poResult.has_orbit
    result.success = false;
    result.message = poResult.message;
    result.msg = poResult.message;
    return
end

orbit = struct();
orbit.t = poResult.orbit_t(:);
orbit.y = poResult.orbit_y;
orbit.period = poResult.period;

observableSpec = @(t, y) y(:, 1) + y(:, 2);
features = evaluate_orbit_features(orbit, observableSpec, [], featureOptions);

result.success = true;
result.message = "";
result.msg = "";
result.orbit = orbit;
result.features = features;
end

function value = get_option_struct(options, name)
value = struct();
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
end
end

function base = merge_struct(base, override)
names = fieldnames(override);
for i = 1:numel(names)
    base.(names{i}) = override.(names{i});
end
end
