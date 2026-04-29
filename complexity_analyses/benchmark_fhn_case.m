function result = benchmark_fhn_case( ...
    case_data, network, base_params, parameter_noise_strength, random_seed, ...
    solver_name, search_window, initial_state_value, ...
    target_property, target_node_indices, target_values, ...
    target_relative_shift, target_absolute_floor, ...
    controlled_parameter_name, controlled_node_indices, ...
    err_bound, max_step_size, continuation_options, newton_options, ...
    assembly_sample_count)
%BENCHMARK_FHN_CASE Run one full-FMAM FHN complexity case.

case_start = tic;
result = initialize_result(case_data, parameter_noise_strength);

validate_modulation_settings( ...
    case_data.N, case_data.targetCount, target_node_indices, target_values, ...
    target_relative_shift, controlled_node_indices);

network.N = case_data.N;
model_seed = random_seed;

model_timer = tic;
model = build_heterogeneous_fhn_network( ...
    case_data.N, network, base_params, parameter_noise_strength, ...
    case_data.netName, case_data.couplingType, initial_state_value, model_seed);
result.modelBuildTime = toc(model_timer);
result.dimVar = model.dimVar;
result.dimParams = model.dimParams;
result.numEdges = nnz(model.omega);
result.networkDensity = nnz(model.omega) / numel(model.omega);

orbit_timer = tic;
orbit_result = networkexp.evaluate_periodic_orbit( ...
    model.rhs, model.initialState, model.rhsParams, solver_name, ...
    search_window, model.solverOptions);
result.orbitTime = toc(orbit_timer);
result.orbitStatus = orbit_result.status;
result.orbitMessage = orbit_result.message;
result.orbitPeriod = orbit_result.observables.period;

if ~orbit_result.success
    result.status = 'orbit_failed';
    result.totalTime = toc(case_start);
    return
end

state_timer = tic;
obs = {};
PV = struct('name', 'var', 'idx', 1);
state_obj = state(obs, model.params, orbit_result.TS{1}, orbit_result.TS{2}, case_data.M, PV);
state_obj.updatePeriod();
state_obj.updateVar2();
solver_view = fmam_state_ops.solverViewFromState(state_obj);
result.stateBuildTime = toc(state_timer);

target_nodes = choose_nodes(target_node_indices, case_data.targetCount, model.N);
controlled_nodes = choose_nodes(controlled_node_indices, numel(target_nodes), model.N);
items_per = build_items_per(state_obj, target_property, target_nodes, target_values, ...
    target_relative_shift, target_absolute_floor);
controlled_idx = model.parameterLayout.(controlled_parameter_name)(controlled_nodes);

constructor_args = {'derivatives', model.derivatives, ...
    'continuationOptions', continuation_options, ...
    'newtonOptions', newton_options, ...
    'verbose', false};
if ~isempty(assembly_sample_count)
    constructor_args = [constructor_args, {'assemblySampleCount', assembly_sample_count}];
end

task = FMAM_ODE(model.system, obs, solver_view, items_per, controlled_idx, ...
    max_step_size, err_bound, constructor_args{:});
task.psiUpdateMode = true;
task.refreshPsiModeReferences();
task.needLog = true;

result.linearSystemSize = estimate_linear_system_size(task.exportSolverView());
result.assemblySampleCount = task.assemblySampleCount;

fit_timer = tic;
fit_result = task.fit();
result.initialFitTime = toc(fit_timer);
result.initialFitConverged = fit_result.converged;
result.initialFitIterations = fit_result.iterations;
result.initialFitResidual = fit_result.scalarError;
if ~fit_result.converged
    result.status = 'initial_fit_failed';
    result.totalTime = toc(case_start);
    return
end

step_timer = tic;
task.step();
result.continuationTime = toc(step_timer);

final_fit_timer = tic;
final_fit_result = task.fit();
result.finalFitTime = toc(final_fit_timer);
result.finalFitConverged = final_fit_result.converged;
result.finalFitIterations = final_fit_result.iterations;
result.finalFitResidual = final_fit_result.scalarError;

result.fmamTime = result.stateBuildTime + result.initialFitTime + ...
    result.continuationTime + result.finalFitTime;
result.acceptedSteps = numel(task.logs);
result.newtonIterationsTotal = sum_numeric_log_field(task.logs, 'newtonIterations') + ...
    result.initialFitIterations + result.finalFitIterations;
result.newtonIterationsPerStepMean = mean_numeric_log_field(task.logs, 'newtonIterations');
result.usedLMCount = sum_logical_log_field(task.logs, 'usedLM');
result.finalLambda = task.continuationStatus.lambda;
result.continuationCompleted = task.continuationStatus.completed;
result.continuationReason = task.continuationStatus.reason;
if ~task.continuationStatus.completed
    result.status = 'continuation_failed';
elseif ~final_fit_result.converged
    result.status = 'final_fit_failed';
else
    result.status = 'success';
end
result.totalTime = toc(case_start);
end

function result = initialize_result(case_data, parameter_noise_strength)
result = struct();
result.caseIndex = case_data.caseIndex;
result.sweepName = case_data.sweepName;
result.N = case_data.N;
result.M = case_data.M;
result.targetCount = case_data.targetCount;
result.netName = case_data.netName;
result.couplingType = case_data.couplingType;
result.errBound = case_data.errBound;
result.dlambdaCap = case_data.dlambdaCap;
result.parameterNoiseStrength = parameter_noise_strength;
result.status = 'not_started';
result.modelBuildTime = NaN;
result.orbitTime = NaN;
result.stateBuildTime = NaN;
result.initialFitTime = NaN;
result.continuationTime = NaN;
result.finalFitTime = NaN;
result.fmamTime = NaN;
result.totalTime = NaN;
result.dimVar = NaN;
result.dimParams = NaN;
result.linearSystemSize = NaN;
result.assemblySampleCount = NaN;
result.numEdges = NaN;
result.networkDensity = NaN;
result.orbitStatus = '';
result.orbitMessage = '';
result.orbitPeriod = NaN;
result.initialFitConverged = false;
result.initialFitIterations = NaN;
result.initialFitResidual = NaN;
result.finalFitConverged = false;
result.finalFitIterations = NaN;
result.finalFitResidual = NaN;
result.acceptedSteps = NaN;
result.newtonIterationsTotal = NaN;
result.newtonIterationsPerStepMean = NaN;
result.usedLMCount = NaN;
result.finalLambda = NaN;
result.continuationCompleted = false;
result.continuationReason = '';
end

function validate_modulation_settings(N, target_count, target_node_indices, target_values, ...
    target_relative_shift, controlled_node_indices)
if target_count > N
    error('fhnComplexity:TooManyTargets', ...
        'targetCount=%d exceeds N=%d. Increase N or reduce target_count_list.', ...
        target_count, N);
end
if ~isempty(target_node_indices) && numel(target_node_indices) < target_count
    error('fhnComplexity:TooFewTargetNodeIndices', ...
        'target_node_indices has %d entries, but targetCount=%d.', ...
        numel(target_node_indices), target_count);
end
if ~isempty(controlled_node_indices) && numel(controlled_node_indices) < target_count
    error('fhnComplexity:TooFewControlledNodeIndices', ...
        'controlled_node_indices has %d entries, but targetCount=%d.', ...
        numel(controlled_node_indices), target_count);
end
if ~isempty(target_values) && numel(target_values) ~= target_count
    error('fhnComplexity:TargetValueCountMismatch', ...
        'target_values must be empty or have targetCount=%d values.', target_count);
end
if isempty(target_values) && ~(isscalar(target_relative_shift) || numel(target_relative_shift) == target_count)
    error('fhnComplexity:TargetShiftCountMismatch', ...
        'target_relative_shift must be scalar or have targetCount=%d values.', target_count);
end
end

function nodes = choose_nodes(node_indices, count, N)
if isempty(node_indices)
    nodes = 1:count;
else
    if numel(node_indices) < count
        error('fhnComplexity:TooFewNodeIndices', ...
            'Only %d node indices were provided, but %d are required.', ...
            numel(node_indices), count);
    end
    nodes = node_indices(1:count);
end
if any(nodes < 1) || any(nodes > N) || any(nodes ~= floor(nodes))
    error('fhnComplexity:InvalidNodeIndices', ...
        'Node indices must be integers between 1 and N=%d.', N);
end
end

function items_per = build_items_per(state_obj, target_property, target_nodes, ...
    target_values, target_relative_shift, target_absolute_floor)
items_per = repmat(struct('prop', target_property, 'idx', 1, 'target', NaN), ...
    1, numel(target_nodes));
base_values = state_obj.(target_property);

if ~isempty(target_values) && numel(target_values) ~= numel(target_nodes)
    error('fhnComplexity:TargetValueCountMismatch', ...
        'target_values must be empty or have one value for each target node.');
end
if isempty(target_values) && ...
        ~(isscalar(target_relative_shift) || numel(target_relative_shift) == numel(target_nodes))
    error('fhnComplexity:TargetShiftCountMismatch', ...
        'target_relative_shift must be scalar or have one value for each target node.');
end

for i = 1:numel(target_nodes)
    idx = target_nodes(i);
    base_value = base_values(idx);
    if isempty(target_values)
        if isscalar(target_relative_shift)
            relative_shift = target_relative_shift;
        else
            relative_shift = target_relative_shift(i);
        end
        shift = relative_shift * max(abs(base_value), target_absolute_floor);
        target_value = base_value + shift;
    else
        target_value = target_values(i);
    end
    items_per(i).prop = target_property;
    items_per(i).idx = idx;
    items_per(i).target = target_value;
end
end

function n = estimate_linear_system_size(solver_view)
n = numel(solver_view.params) + numel(solver_view.p_Psi) + ...
    numel(solver_view.q_Psi) + numel(solver_view.p_var) + ...
    numel(solver_view.q_var) + numel(solver_view.varPhiMax) + ...
    numel(solver_view.varPhiMin) + numel(solver_view.obsPhiMax) + ...
    numel(solver_view.obsPhiMin);
end

function value = sum_numeric_log_field(logs, field_name)
if isempty(logs)
    value = 0;
    return
end
values = [logs.(field_name)];
value = sum(values(isfinite(values)));
end

function value = mean_numeric_log_field(logs, field_name)
if isempty(logs)
    value = NaN;
    return
end
values = [logs.(field_name)];
values = values(isfinite(values));
value = mean(values);
end

function value = sum_logical_log_field(logs, field_name)
if isempty(logs)
    value = 0;
    return
end
value = sum(logical([logs.(field_name)]));
end
