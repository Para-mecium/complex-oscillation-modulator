function model = build_model(cfg)
circadian.ensure_paths();

paramNames = {'Kd', 'AT'};
constants = cfg.clock.constants;
defaultParams = cfg.clock.initialParams;
system = build_clock_system(constants);
observables = {@(variable) variable(:, 2) + variable(:, 3)};
pv = struct('name', 'obs', 'idx', 1);
% pv = struct('name', 'var', 'idx', 2);
derivatives = get_derivative_cache(constants, system, observables, numel(paramNames));

model = struct();
model.kind = 'circadian';
model.paramNames = paramNames;
model.constants = constants;
model.defaultParams = reshape(defaultParams, 1, []);
model.system = system;
model.observables = observables;
model.derivatives = derivatives;
model.pv = pv;
model.rhs = @(y, p) rhs_vector(system, y, p);
end

function system = build_clock_system(c)
system = { ...
    @(variable, parameter) c.beta * (free_active(parameter(2), parameter(1), variable(:, 3)) ./ parameter(2) - variable(:, 1)), ...
    @(variable, parameter) c.beta * (variable(:, 1) - variable(:, 2)), ...
    @(variable, parameter) c.beta * (variable(:, 2) - variable(:, 3))};
end

function value = free_active(AT, Kd, Pn)
delta = AT - Pn - Kd;
value = 0.5 * (delta + sqrt(delta.^2 + 4 * Kd * AT));
end

function derivatives = get_derivative_cache(constants, system, observables, nParams)
persistent cache
if isempty(cache)
    cache = struct();
end

key = matlab.lang.makeValidName(signature_from_constants(constants));
if isfield(cache, key)
    derivatives = cache.(key);
    return
end

derivatives = build_symbolic_derivatives(system, observables, nParams);
cache.(key) = derivatives;
end

function signature = signature_from_constants(constants)
names = sort(fieldnames(constants));
parts = cell(1, numel(names));
for i = 1:numel(names)
    parts{i} = sprintf('%s_%s', names{i}, mat2str(constants.(names{i}), 12));
end
signature = strjoin(parts, '__');
end

function rhs = rhs_vector(system, y, p)
stateRow = reshape(y, 1, []);
paramRow = reshape(p, 1, []);
rhs = zeros(numel(system), 1);
for i = 1:numel(system)
    rhs(i) = system{i}(stateRow, paramRow);
end
end
