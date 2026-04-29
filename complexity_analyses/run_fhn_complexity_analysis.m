clear
clc

%% Paths
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(script_dir);
addpath(repo_dir);
addpath(script_dir);
addpath(fullfile(repo_dir, 'network_modulatability'));
addpath(fullfile(repo_dir, 'PO_extract'));

%% User settings
random_seed = 1;
output_dir = fullfile(script_dir, 'data');
save_intermediate = true;

% Complexity variables.
N_list = [5 10 20 30 40 50 75 100];
M_list = [5 10 20 30 40 50];
target_count_list = [1 2 5 10 15 20];
err_bound_list = [1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10];
dlambda_cap_list = [0.005 0.01 0.02 0.05 0.1];
network_type_list = {'FC'};
coupling_type_list = {'synapse'};

% Baseline case used when sweeping one complexity variable at a time.
baseline_N = 20;
baseline_M = 10;
baseline_target_count = 5;
baseline_network_type = 'FC';
baseline_coupling_type = 'synapse';
baseline_err_bound = 1e-6;
baseline_dlambda_cap = 0.05;

run_dimension_sweep = true;
run_truncation_sweep = true;
run_target_count_sweep = true;
run_err_bound_sweep = true;
run_dlambda_cap_sweep = true;
run_network_coupling_sweep = false;

% Network definition. ER/BA/SW follow networkexp.build_network_matrix.
% FC is added here: all off-diagonal edges are present, with U[0.2, 1] weights.
network = struct();
network.degAvg = 5;
network.K = 0.1;
network.baSeedDegree = 5;
network.swNeighborCount = 2;
network.swRewireProb = 0.15;

% FHN parameters. Scalars are expanded to length N; vectors can be supplied
% later for neuron-specific tuning.
base_params = struct();
base_params.theta = 0.5;
base_params.gamma = 1;
base_params.epsilon = 5e-2;
base_params.tau = 1;
base_params.I = 0.38;
base_params.V = 0;
parameter_noise_strength = 0;

% Coupling types:
% synapse: omega * sigma, sigma = 1 ./ (1 + exp(-tau .* (v - V)))
% gap:     omega * v - sum(omega, 2) .* v

% FMAM target and controlled parameter.
target_property = 'varMax';
target_node_indices = [];
target_values = [];
target_relative_shift = 0.01;
target_absolute_floor = 1e-3;

controlled_parameter_name = 'I';
controlled_node_indices = [];

% Numerical settings.
solver_name = 'ode15s';
search_window = 5000;
initial_state_value = 0.21;
max_step_size = [];
assembly_sample_count = [];

continuation_options = struct( ...
    'initialLambdaStep', baseline_dlambda_cap, ...
    'maxLambdaStep', baseline_dlambda_cap, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);
newton_options = struct();

%% Build benchmark cases
case_count = 0;
if run_dimension_sweep
    case_count = case_count + numel(N_list);
end
if run_truncation_sweep
    case_count = case_count + numel(M_list);
end
if run_target_count_sweep
    case_count = case_count + numel(target_count_list);
end
if run_err_bound_sweep
    case_count = case_count + numel(err_bound_list);
end
if run_dlambda_cap_sweep
    case_count = case_count + numel(dlambda_cap_list);
end
if run_network_coupling_sweep
    case_count = case_count + numel(network_type_list) * numel(coupling_type_list);
end
if case_count == 0
    error('fhnComplexity:NoCasesSelected', ...
        'At least one sweep flag must be true.');
end

case_rows = repmat(make_case('', NaN, NaN, NaN, '', '', NaN, NaN), 1, case_count);
case_cursor = 0;

if run_dimension_sweep
    for N = N_list
        case_cursor = case_cursor + 1;
        case_rows(case_cursor) = make_case('dimension', N, baseline_M, ...
            baseline_target_count, baseline_network_type, baseline_coupling_type, ...
            baseline_err_bound, baseline_dlambda_cap);
    end
end

if run_truncation_sweep
    for M = M_list
        case_cursor = case_cursor + 1;
        case_rows(case_cursor) = make_case('truncation', baseline_N, M, ...
            baseline_target_count, baseline_network_type, baseline_coupling_type, ...
            baseline_err_bound, baseline_dlambda_cap);
    end
end

if run_target_count_sweep
    for target_count = target_count_list
        case_cursor = case_cursor + 1;
        case_rows(case_cursor) = make_case('targets', baseline_N, baseline_M, ...
            target_count, baseline_network_type, baseline_coupling_type, ...
            baseline_err_bound, baseline_dlambda_cap);
    end
end

if run_err_bound_sweep
    for err_bound = err_bound_list
        case_cursor = case_cursor + 1;
        case_rows(case_cursor) = make_case('errBound', baseline_N, baseline_M, ...
            baseline_target_count, baseline_network_type, baseline_coupling_type, ...
            err_bound, baseline_dlambda_cap);
    end
end

if run_dlambda_cap_sweep
    for dlambda_cap = dlambda_cap_list
        case_cursor = case_cursor + 1;
        case_rows(case_cursor) = make_case('dlambdaCap', baseline_N, baseline_M, ...
            baseline_target_count, baseline_network_type, baseline_coupling_type, ...
            baseline_err_bound, dlambda_cap);
    end
end

if run_network_coupling_sweep
    for i_net = 1:numel(network_type_list)
        for i_coupling = 1:numel(coupling_type_list)
            case_cursor = case_cursor + 1;
            case_rows(case_cursor) = make_case('network_coupling', baseline_N, baseline_M, ...
                baseline_target_count, network_type_list{i_net}, coupling_type_list{i_coupling}, ...
                baseline_err_bound, baseline_dlambda_cap);
        end
    end
end

for i_case = 1:numel(case_rows)
    case_rows(i_case).caseIndex = i_case;
end

%% Run benchmarks
baseline_dir = build_baseline_dir_label( ...
    baseline_N, baseline_M);

for i_case = 1:numel(case_rows)
    case_data = case_rows(i_case);
    fprintf(['[%d/%d] %s: N=%d, M=%d, targets=%d, net=%s, coupling=%s, ' ...
        'errBound=%g, dlambdaCap=%g\n'], ...
        i_case, numel(case_rows), case_data.sweepName, case_data.N, case_data.M, ...
        case_data.targetCount, case_data.netName, case_data.couplingType, ...
        case_data.errBound, case_data.dlambdaCap);

    case_continuation_options = continuation_options;
    case_continuation_options.initialLambdaStep = case_data.dlambdaCap;
    case_continuation_options.maxLambdaStep = case_data.dlambdaCap;

    result = benchmark_fhn_case( ...
        case_data, network, base_params, parameter_noise_strength, random_seed, ...
        solver_name, search_window, initial_state_value, ...
        target_property, target_node_indices, target_values, ...
        target_relative_shift, target_absolute_floor, ...
        controlled_parameter_name, controlled_node_indices, ...
        case_data.errBound, max_step_size, case_continuation_options, newton_options, ...
        assembly_sample_count);

    fprintf('  status=%s, reason=%s, lambda=%.4g, FMAM %.3fs, steps=%g, Newton=%g\n', ...
        result.status, result.continuationReason, ...
        result.finalLambda, result.fmamTime, ...
        result.acceptedSteps, result.newtonIterationsTotal);

    [case_output_dir, mat_file] = build_case_output_path(output_dir, baseline_dir, case_data);
    if exist(case_output_dir, 'dir') ~= 7
        mkdir(case_output_dir);
    end
    baseline_settings = struct( ...
        'N', baseline_N, ...
        'M', baseline_M, ...
        'targetCount', baseline_target_count, ...
        'netName', baseline_network_type, ...
        'couplingType', baseline_coupling_type, ...
        'errBound', baseline_err_bound, ...
        'dlambdaCap', baseline_dlambda_cap, ...
        'targetProperty', target_property, ...
        'controlledParameter', controlled_parameter_name, ...
        'randomSeed', random_seed, ...
        'networkK', network.K);
    save(mat_file, 'case_data', 'result', 'baseline_settings', ...
        'network', 'base_params', 'parameter_noise_strength', ...
        'target_property', 'target_node_indices', 'target_values', ...
        'target_relative_shift', 'target_absolute_floor', ...
        'controlled_parameter_name', 'controlled_node_indices', ...
        'solver_name', 'search_window', 'initial_state_value', ...
        'max_step_size', 'assembly_sample_count', ...
        'case_continuation_options', 'newton_options');
    fprintf('  saved: %s\n', mat_file);
end

fprintf('Saved case MAT files under: %s\n', fullfile(output_dir, baseline_dir));

function row = make_case(sweep_name, N, M, target_count, net_name, coupling_type, ...
    err_bound, dlambda_cap)
row = struct( ...
    'caseIndex', NaN, ...
    'sweepName', sweep_name, ...
    'N', N, ...
    'M', M, ...
    'targetCount', target_count, ...
    'netName', net_name, ...
    'couplingType', coupling_type, ...
    'errBound', err_bound, ...
    'dlambdaCap', dlambda_cap);
end

function label = build_baseline_dir_label(baseline_N, baseline_M)
label = sprintf('base_N%d_M%d', baseline_N, baseline_M);
label = sanitize_path_part(label);
end

function [case_output_dir, mat_file] = build_case_output_path(output_dir, baseline_dir, case_data)
switch case_data.sweepName
    case 'dimension'
        sweep_dir = 'sys_dim_sweep';
        file_name = sprintf('sys_dim=%s.mat', value_label(case_data.N));
    case 'truncation'
        sweep_dir = 'truncation_order_sweep';
        file_name = sprintf('truncation_order=%s.mat', value_label(case_data.M));
    case 'targets'
        sweep_dir = 'target_num_sweep';
        file_name = sprintf('target_num=%s.mat', value_label(case_data.targetCount));
    case 'errBound'
        sweep_dir = 'err_bound_sweep';
        file_name = sprintf('err_bound=%s.mat', value_label(case_data.errBound));
    case 'dlambdaCap'
        sweep_dir = 'dlambda_cap_sweep';
        file_name = sprintf('dlambda_cap=%s.mat', value_label(case_data.dlambdaCap));
    case 'network_coupling'
        sweep_dir = 'network_coupling_sweep';
        file_name = sprintf('network=%s_coupling=%s.mat', ...
            case_data.netName, case_data.couplingType);
    otherwise
        error('fhnComplexity:UnknownSweepName', ...
            'Unsupported sweep name: %s.', case_data.sweepName);
end
case_output_dir = fullfile(output_dir, baseline_dir, sweep_dir);
mat_file = fullfile(case_output_dir, sanitize_file_name(file_name));
end

function text = value_label(value)
text = sprintf('%g', value);
text = strrep(text, '.', 'p');
text = strrep(text, '-', 'm');
end

function text = sanitize_path_part(text)
text = regexprep(text, '[^A-Za-z0-9_-]', '');
end

function text = sanitize_file_name(text)
text = regexprep(text, '[^A-Za-z0-9_=.-]', '');
end
