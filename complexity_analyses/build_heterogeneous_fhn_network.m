function model = build_heterogeneous_fhn_network( ...
    N, network, base_params, parameter_noise_strength, net_name, coupling_type, ...
    initial_state_value, seed)
%BUILD_HETEROGENEOUS_FHN_NETWORK Build vector-parameter FHN network for FMAM.

network.N = N;
if strcmpi(net_name, 'FC')
    [omega, network_metadata] = build_fully_connected_network(network, seed);
else
    [omega, network_metadata] = networkexp.build_network_matrix(network, net_name, seed);
end

layout = make_parameter_layout(N);
params = make_parameter_vector(base_params, parameter_noise_strength, N, layout, seed);
spec = struct('N', N, 'omega', sparse(omega), ...
    'couplingType', coupling_type, 'layout', layout);

model = struct();
model.N = N;
model.dimVar = 2 * N;
model.dimParams = numel(params);
model.params = params;
model.parameterLayout = layout;
model.omega = sparse(omega);
model.networkMetadata = network_metadata;
model.couplingType = coupling_type;
model.initialState = initial_state_value * ones(2 * N, 1);
model.rhsParams = struct('params', params, 'spec', spec);
model.rhs = @fhn_rhs_for_orbit;
model.system = build_fmam_system(spec);
model.derivatives = build_fmam_derivatives(spec);
model.solverOptions = build_solver_options(spec, params);
end

function layout = make_parameter_layout(N)
layout.theta = 1:N;
layout.gamma = N + (1:N);
layout.epsilon = 2 * N + (1:N);
layout.tau = 3 * N + (1:N);
layout.I = 4 * N + (1:N);
layout.V = 5 * N + (1:N);
layout.names = {'theta', 'gamma', 'epsilon', 'tau', 'I', 'V'};
end

function [omega, metadata] = build_fully_connected_network(network, seed)
previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng));
rng(seed, 'twister');

N = network.N;
weights = 0.2 + 0.8 * rand(N, N);
weights(1:N + 1:end) = 0;
omega = sparse((network.K / N) * weights);

metadata = struct( ...
    'netName', 'FC', ...
    'seed', seed, ...
    'numNodes', N, ...
    'degreeAverage', N - 1, ...
    'couplingScale', network.K, ...
    'numEdges', nnz(omega), ...
    'density', nnz(omega) / numel(omega));
end

function params = make_parameter_vector(base_params, noise_strength, N, layout, seed)
previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng));
rng(seed + 7919, 'twister');

params = zeros(1, 6 * N);
for i = 1:numel(layout.names)
    name = layout.names{i};
    values = base_params.(name);
    if isscalar(values)
        values = values * ones(1, N);
    else
        values = reshape(values, 1, []);
        if numel(values) ~= N
            error('fhnComplexity:ParameterSizeMismatch', ...
                'base_params.%s must be scalar or length N=%d.', name, N);
        end
    end

    if noise_strength ~= 0
        noise = 2 * rand(1, N) - 1;
        if strcmp(name, 'V')
            values = values + noise_strength * noise;
        else
            values = values .* (1 + noise_strength * noise);
        end
    end
    params(layout.(name)) = values;
end
end

function system = build_fmam_system(spec)
system = cell(1, 2 * spec.N);
for eq_idx = 1:(2 * spec.N)
    system{eq_idx} = @(variable, params) fhn_system_component(variable, params, eq_idx, spec);
end
end

function derivatives = build_fmam_derivatives(spec)
dim_var = 2 * spec.N;
dim_params = 6 * spec.N;
zero_derivative = @(variable, params) zeros(size(variable, 1), 1);
derivatives = struct();
derivatives.var = repmat(struct('function', zero_derivative), dim_var, dim_var + dim_params);
derivatives.obs = [];
derivatives.obs2 = [];

N = spec.N;
omega_pattern = spones(spec.omega);
for i = 1:N
    v_eq = i;
    w_eq = N + i;

    variable_idx = unique([i, N + i, find(omega_pattern(i, :))]);
    param_idx = unique([spec.layout.theta(i), spec.layout.gamma(i), ...
        spec.layout.epsilon(i), spec.layout.I(i)]);
    if strcmp(spec.couplingType, 'synapse')
        source_nodes = find(omega_pattern(i, :));
        param_idx = unique([param_idx, spec.layout.tau(source_nodes), spec.layout.V(source_nodes)]);
    end

    for idx = variable_idx
        derivatives.var(v_eq, idx).function = ...
            make_derivative_handle(v_eq, idx, spec);
    end
    for idx = param_idx
        derivatives.var(v_eq, dim_var + idx).function = ...
            make_derivative_handle(v_eq, dim_var + idx, spec);
    end

    for idx = [i, N + i]
        derivatives.var(w_eq, idx).function = ...
            make_derivative_handle(w_eq, idx, spec);
    end
    for idx = [spec.layout.gamma(i), spec.layout.epsilon(i)]
        derivatives.var(w_eq, dim_var + idx).function = ...
            make_derivative_handle(w_eq, dim_var + idx, spec);
    end
end
end

function handle = make_derivative_handle(eq_idx, diff_idx, spec)
handle = @(variable, params) fhn_derivative_component(variable, params, eq_idx, diff_idx, spec);
end

function solver_options = build_solver_options(spec, params)
N = spec.N;
identity = speye(N);
omega_pattern = spones(spec.omega);
solver_options = struct();
solver_options.JPattern = [omega_pattern + identity, identity; identity, identity];
solver_options.Jacobian = @(t, x) fhn_state_jacobian(x, params, spec);
end

function dydt = fhn_rhs_for_orbit(~, x, rhs_params)
spec = rhs_params.spec;
params = rhs_params.params;
dydt = zeros(2 * spec.N, 1);
variable = reshape(x, 1, []);
for eq_idx = 1:(2 * spec.N)
    dydt(eq_idx) = fhn_system_component(variable, params, eq_idx, spec);
end
end

function value = fhn_system_component(variable, params, eq_idx, spec)
N = spec.N;
theta = params(spec.layout.theta);
gamma = params(spec.layout.gamma);
epsilon = params(spec.layout.epsilon);
tau = params(spec.layout.tau);
current_I = params(spec.layout.I);
threshold_V = params(spec.layout.V);

v = variable(:, 1:N);
w = variable(:, N + (1:N));

if eq_idx <= N
    i = eq_idx;
    local = v(:, i) .* (v(:, i) - theta(i)) .* (1 - v(:, i)) - w(:, i) + current_I(i);
    value = local + coupling_component(v, tau, threshold_V, i, spec);
else
    i = eq_idx - N;
    value = epsilon(i) .* (v(:, i) - gamma(i) .* w(:, i));
end
end

function value = coupling_component(v, tau, threshold_V, i, spec)
row = spec.omega(i, :);
if strcmp(spec.couplingType, 'synapse')
    sigma = 1 ./ (1 + exp(-tau .* (v - threshold_V)));
    value = full(sigma * row');
elseif strcmp(spec.couplingType, 'gap')
    value = full(v * row') - full(sum(row, 2)) .* v(:, i);
else
    error('fhnComplexity:UnknownCouplingType', ...
        'coupling_type must be ''synapse'' or ''gap''.');
end
end

function value = fhn_derivative_component(variable, params, eq_idx, diff_idx, spec)
N = spec.N;
dim_var = 2 * N;
is_param = diff_idx > dim_var;
param_idx = diff_idx - dim_var;

theta = params(spec.layout.theta);
gamma = params(spec.layout.gamma);
epsilon = params(spec.layout.epsilon);
tau = params(spec.layout.tau);
threshold_V = params(spec.layout.V);

v = variable(:, 1:N);
w = variable(:, N + (1:N));
value = zeros(size(variable, 1), 1);

if eq_idx <= N
    i = eq_idx;
    if ~is_param
        if diff_idx <= N
            j = diff_idx;
            if j == i
                value = value + (-3 * v(:, i).^2 + 2 * (1 + theta(i)) .* v(:, i) - theta(i));
            end
            value = value + coupling_derivative_by_voltage(variable, params, i, j, spec);
        elseif diff_idx == N + i
            value = value - 1;
        end
        return
    end

    if param_idx == spec.layout.theta(i)
        value = -v(:, i) .* (1 - v(:, i));
    elseif param_idx == spec.layout.I(i)
        value = ones(size(variable, 1), 1);
    elseif strcmp(spec.couplingType, 'synapse')
        source_node = find_synapse_parameter_node(param_idx, spec.layout.tau, spec.layout.V);
        if source_node > 0
            value = synapse_parameter_derivative(v, tau, threshold_V, i, source_node, param_idx, spec);
        end
    end
else
    i = eq_idx - N;
    if ~is_param
        if diff_idx == i
            value = epsilon(i) * ones(size(variable, 1), 1);
        elseif diff_idx == N + i
            value = -epsilon(i) * gamma(i) * ones(size(variable, 1), 1);
        end
        return
    end

    if param_idx == spec.layout.gamma(i)
        value = -epsilon(i) .* w(:, i);
    elseif param_idx == spec.layout.epsilon(i)
        value = v(:, i) - gamma(i) .* w(:, i);
    end
end
end

function value = coupling_derivative_by_voltage(variable, params, target_node, source_node, spec)
N = spec.N;
v = variable(:, 1:N);
tau = params(spec.layout.tau);
threshold_V = params(spec.layout.V);
weight = full(spec.omega(target_node, source_node));

if strcmp(spec.couplingType, 'synapse')
    sigma = 1 ./ (1 + exp(-tau(source_node) .* (v(:, source_node) - threshold_V(source_node))));
    value = weight .* tau(source_node) .* sigma .* (1 - sigma);
elseif strcmp(spec.couplingType, 'gap')
    value = weight * ones(size(variable, 1), 1);
    if source_node == target_node
        value = value - full(sum(spec.omega(target_node, :), 2));
    end
else
    error('fhnComplexity:UnknownCouplingType', ...
        'coupling_type must be ''synapse'' or ''gap''.');
end
end

function node = find_synapse_parameter_node(param_idx, tau_layout, V_layout)
node = find(tau_layout == param_idx, 1);
if isempty(node)
    node = find(V_layout == param_idx, 1);
end
if isempty(node)
    node = 0;
end
end

function value = synapse_parameter_derivative(v, tau, threshold_V, ...
    target_node, source_node, param_idx, spec)
weight = full(spec.omega(target_node, source_node));
sigma = 1 ./ (1 + exp(-tau(source_node) .* (v(:, source_node) - threshold_V(source_node))));
sigma_prime = sigma .* (1 - sigma);

if param_idx == spec.layout.tau(source_node)
    value = weight .* (v(:, source_node) - threshold_V(source_node)) .* sigma_prime;
else
    value = -weight .* tau(source_node) .* sigma_prime;
end
end

function jacobian = fhn_state_jacobian(x, params, spec)
N = spec.N;
variable = reshape(x, 1, []);
max_entries = (2 * N)^2;
rows = zeros(max_entries, 1);
cols = zeros(max_entries, 1);
vals = zeros(max_entries, 1);
entry_count = 0;
for eq_idx = 1:(2 * N)
    for var_idx = 1:(2 * N)
        val = fhn_derivative_component(variable, params, eq_idx, var_idx, spec);
        if val ~= 0
            entry_count = entry_count + 1;
            rows(entry_count) = eq_idx;
            cols(entry_count) = var_idx;
            vals(entry_count) = val;
        end
    end
end
jacobian = sparse(rows(1:entry_count), cols(1:entry_count), vals(1:entry_count), 2 * N, 2 * N);
end
