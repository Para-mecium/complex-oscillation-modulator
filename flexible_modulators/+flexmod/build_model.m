function model = build_model(cfg)
flexmod.ensure_paths();

switch lower(cfg.modelType)
    case 'base'
        paramNames = {'I1', 'ET'};
        constants = cfg.base.constants;
        defaultParams = cfg.base.initialParams;
        system = build_base_system(constants);
    case 'temperature'
        paramNames = {'ET', 'Temperature'};
        constants = cfg.temperature.constants;
        defaultParams = cfg.temperature.initialParams;
        system = build_temperature_system(constants);
    otherwise
        error('flexmod:UnknownModelType', 'Unsupported modelType ''%s''.', cfg.modelType);
end

observables = {};
pv = struct('name', 'var', 'idx', 2);
derivatives = get_derivative_cache(cfg.modelType, constants, system, observables, numel(paramNames));

model = struct();
model.kind = lower(cfg.modelType);
model.paramNames = paramNames;
model.constants = constants;
model.defaultParams = reshape(defaultParams, 1, []);
model.system = system;
model.observables = observables;
model.derivatives = derivatives;
model.pv = pv;
model.rhs = @(y, p) rhs_vector(system, y, p);
end

function system = build_base_system(c)
system = { ...
    @(variable, parameter) c.k1 * c.S * c.Kd^c.p ./ (c.Kd^c.p + variable(:,2).^c.p) * ...
        regulated_drive(parameter(1), c) - c.kdx * variable(:,1), ...
    @(variable, parameter) c.ksy * variable(:,1) - c.kdy * variable(:,2) - ...
        c.k2 * parameter(2) * variable(:,2) ./ (c.Km + variable(:,2) + c.KI * variable(:,2).^2)};
end

function system = build_temperature_system(c)
system = { ...
    @(variable, parameter) arrhenius_rate(c.k1_ref, c.Ek1, parameter(2), c) * ...
        c.S * c.Kd^c.p ./ (c.Kd^c.p + variable(:,2).^c.p) .* fixed_drive(c) - ...
        arrhenius_rate(c.kdx_ref, c.Ekdx, parameter(2), c) * variable(:,1), ...
    @(variable, parameter) arrhenius_rate(c.ksy_ref, c.Eksy, parameter(2), c) * variable(:,1) - ...
        arrhenius_rate(c.kdy_ref, c.Ekdy, parameter(2), c) * variable(:,2) - ...
        arrhenius_rate(c.k2_ref, c.Ek2, parameter(2), c) * parameter(1) * ...
        variable(:,2) ./ (c.Km + variable(:,2) + c.KI * variable(:,2).^2)};
end

function value = regulated_drive(I1, c)
hU = c.U * I1^2 / (c.K1 + I1^2);
value = c.bU * hU^2 / (c.KU + hU^2);
end

function value = fixed_drive(c)
hU = c.U * c.I1^2 / (c.K1 + c.I1^2);
value = c.bU * hU^2 / (c.KU + hU^2);
end

function rate = arrhenius_rate(baseRate, activationEnergy, temperature, c)
rate = baseRate * exp(activationEnergy / c.R * (1 / c.Tref - 1 / temperature));
end

function derivatives = get_derivative_cache(modelType, constants, system, observables, nParams)
persistent cache
if isempty(cache)
    cache = struct();
end

key = matlab.lang.makeValidName([modelType '_' signature_from_constants(constants)]);
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
    value = constants.(names{i});
    parts{i} = sprintf('%s_%s', names{i}, mat2str(value, 12));
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
