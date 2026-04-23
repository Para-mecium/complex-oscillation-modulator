clear
% clc

%% Local configuration
need_path = false;
N = 2;
M = 75;

% Fixed network coupling.
seed = 1;
rng(seed, 'twister');

couplingType = 'gap';
G_scale = 0.05;
G = G_scale * (2*randn(N, N)-1);

% Fixed HH parameters.
C0 = 1.0 * ones(N, 1);
gNa0 = 120.0 * ones(N, 1);
ENa0 = 50.0 * ones(N, 1);
gK0 = 36.0 * ones(N, 1);
EK0 = -77.0 * ones(N, 1);
gL0 = 0.3 * ones(N, 1);
EL0 = -54.387 * ones(N, 1);
tau0 = 1.0 * ones(N, 1);
Vstar0 = 0 * ones(N, 1);
Esyn0 = 0.0 * ones(N, 1);
I0_vector = 120 * ones(N, 1);

% Active parameters selected by name.
active_param_names = {'I_1','I_2'};

% Shared initial state.
V0 = 0.2;
m0 = 0.5;
h0 = 0.2;
n0 = 0.5;

% Target amplitudes relative to the initial periodic orbit.
targetScale1 = 4;
targetScale2 = 4;

% Orbit extraction settings.
orbitOptions = struct();
orbitOptions.poOptions = struct( ...
    'solver_name', 'ode15s', ...
    'single_timespan', 1200, ...
    'max_windows', 5, ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 6, ...
    'transientFraction', 0.5, ...
    'samplesPerCycle', 500, ...
    'extractNumPoints', 1200);

% FMAM continuation settings.
errBound = 1e-6;
tightErrBound = 1e-10;
if need_path
    initialLambdaStep = 0.01;
else
    initialLambdaStep = 0.1;
end
continuationOptions = struct( ...
    'initialLambdaStep', initialLambdaStep, ...
    'predictorMode', 'constant',...
    'minLambdaStep', 1e-6);
newtonOptions = struct(...
    'directConditionThreshold', 1e-12,...
    'lmConditionThreshold', 1e-14);

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
addpath(repoDir, '-begin');
addpath(scriptDir, '-begin');
addpath(fullfile(repoDir, 'PO_extract'), '-begin');

%% Build fixed HH parameter struct
defaultParams = struct();
defaultParams.couplingType = couplingType;
defaultParams.I = I0_vector(:);
defaultParams.C = C0(:);
defaultParams.gNa = gNa0(:);
defaultParams.ENa = ENa0(:);
defaultParams.gK = gK0(:);
defaultParams.EK = EK0(:);
defaultParams.gL = gL0(:);
defaultParams.EL = EL0(:);
defaultParams.tau = tau0(:);
defaultParams.Vstar = Vstar0(:);
defaultParams.Esyn = Esyn0(:);
defaultParams.G = G;

y0 = [ ...
    V0 * randn(N, 1); ...
    m0 * ones(N, 1); ...
    h0 * ones(N, 1); ...
    n0 * ones(N, 1)];

%% Initial periodic orbit and features
initialOrbitResult = HH_forward_orbit(build_hh_p(defaultParams), defaultParams.I, y0, N, orbitOptions);

active_specs = build_active_specs(active_param_names);
params = build_active_param_vector(defaultParams, active_specs);
t = initialOrbitResult.orbit.t;
TS_var = initialOrbitResult.orbit.y;
obs = [];

%% Build HH system for FMAM_ODE
sys = build_hh_system(defaultParams, active_specs, N);
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

PV = struct('name', 'var', 'idx', 1);
discretization = fmam_state_defaults.defaultDiscretization();
discretization.reconstruction.phaseMode = 'linearTime';
State = state(obs, params, t, TS_var, M, PV, 'discretization', discretization);
StateView = fmam_state_ops.solverViewFromState(State);

initialfeat = reshape(initialOrbitResult.features.state.amplitude, 1, []);
targetfeat = [ ...
    targetScale1 * initialfeat(1), ...
    targetScale2 * initialfeat(2)];
% targetfeat = [ ...
%     targetScale1 * initialfeat(1)];

items_per = struct;
items_per(1).prop = 'varMax';
items_per(1).idx = 1;
items_per(1).target = 40;

items_per(2).prop = 'varMax';
items_per(2).idx = 2;
items_per(2).target = 32;

% items_per(3).prop = 'varAmp';
% items_per(3).idx = 3;
% items_per(3).target = targetfeat(3);

items_controlled = 1:numel(params);
task = FMAM_ODE(sys, obs, StateView, items_per, items_controlled, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions,'newtonOptions',newtonOptions);
task.psiUpdateMode = false;
task.refreshPsiModeReferences();
task.needLog = need_path;

%% Run modulation
tic
fitResult1 = task.fit();
task.step()
% task.errBound = tightErrBound;
% fitResult2 = task.fit();
% task.step()
elapsedTime = toc;
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

%% Export modulated state and re-extract the ODE periodic orbit
StateView_modulated = task.exportSolverView();
StateDerived_modulated = task.exportDerivedView();

modulatedParams = apply_active_params(defaultParams, active_specs, StateView_modulated.params(:).');
I_modulated = modulatedParams.I(:);
y0_modulated = StateDerived_modulated.TS_var(1, :).';

finalOrbitOptions = struct();
finalOrbitOptions.poOptions = struct( ...
    'solver_name', 'ode15s', ...
    'tspan', [0, max(100 * StateDerived_modulated.period, StateDerived_modulated.t(end))], ...
    'event', 2, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0, ...
    'samplesPerCycle', 500, ...
    'extractNumPoints', 600);

finalOrbitResult = HH_forward_orbit(build_hh_p(modulatedParams), I_modulated, y0_modulated, N, finalOrbitOptions);
finalVarAmp = reshape(finalOrbitResult.features.state.amplitude, 1, []);

%% Report
initialTargetValues = extract_target_feature_values(items_per, initialOrbitResult.features, params);
modulatedTargetValues = extract_target_feature_values(items_per, finalOrbitResult.features, StateView_modulated.params(:).');

disp('Initial active parameters')
for k = 1:numel(active_param_names)
    fprintf('%s=%g\n', active_param_names{k}, params(k));
end

disp('Modulated active parameters')
for k = 1:numel(active_param_names)
    fprintf('%s=%g\n', active_param_names{k}, StateView_modulated.params(k));
end

disp('Initial features')
for k = 1:numel(items_per)
    fprintf('%s=%g\n', target_feature_name(items_per(k)), initialTargetValues(k));
end

disp('Modulated features')
for k = 1:numel(items_per)
    fprintf('%s=%g\n', target_feature_name(items_per(k)), modulatedTargetValues(k));
end

%% Visualization
figure
subplot(2, 1, 1)
plot(initialOrbitResult.orbit.t, initialOrbitResult.orbit.y(:, 1), 'LineWidth', 1.2)
hold on
plot(initialOrbitResult.orbit.t, initialOrbitResult.orbit.y(:, 2), 'LineWidth', 1.2)
hold off
grid on
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Initial HH Periodic Orbit')
legend('V_1', 'V_2', 'Location', 'best')

subplot(2, 1, 2)
plot(finalOrbitResult.orbit.t, finalOrbitResult.orbit.y(:, 1), 'LineWidth', 1.2)
hold on
plot(finalOrbitResult.orbit.t, finalOrbitResult.orbit.y(:, 2), 'LineWidth', 1.2)
hold off
grid on
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Modulated HH Periodic Orbit')
legend('V_1', 'V_2', 'Location', 'best')

function sys = build_hh_system(defaultParams, active_specs, N)
sys = cell(1, 4 * N);
for i = 1:N
    sys{i} = @(variable, parameter) hh_voltage_rhs(variable, parameter, defaultParams, active_specs, N, i);
    sys{N + i} = @(variable, parameter) hh_m_rhs(variable, parameter, defaultParams, active_specs, N, i);
    sys{2 * N + i} = @(variable, parameter) hh_h_rhs(variable, parameter, defaultParams, active_specs, N, i);
    sys{3 * N + i} = @(variable, parameter) hh_n_rhs(variable, parameter, defaultParams, active_specs, N, i);
end
end

function value = hh_voltage_rhs(variable, parameter, defaultParams, active_specs, N, i)
V = variable(:,1:N);
m = variable(:,N + 1:2 * N);
h = variable(:,2 * N + 1:3 * N);
n = variable(:,3 * N + 1:4 * N);

I_i = get_node_param('I', i, parameter, defaultParams, active_specs);
C_i = get_node_param('C', i, parameter, defaultParams, active_specs);
gNa_i = get_node_param('gNa', i, parameter, defaultParams, active_specs);
ENa_i = get_node_param('ENa', i, parameter, defaultParams, active_specs);
gK_i = get_node_param('gK', i, parameter, defaultParams, active_specs);
EK_i = get_node_param('EK', i, parameter, defaultParams, active_specs);
gL_i = get_node_param('gL', i, parameter, defaultParams, active_specs);
EL_i = get_node_param('EL', i, parameter, defaultParams, active_specs);
Esyn_i = get_node_param('Esyn', i, parameter, defaultParams, active_specs);

I_cpl = 0;
switch lower(string(defaultParams.couplingType))
    case "gap"
        for j = 1:N
            if j ~= i
                G_ij = get_matrix_param('G', i, j, parameter, defaultParams, active_specs);
                I_cpl = I_cpl + G_ij .* (V(:,j) - V(:,i));
            end
        end
    case "synapse"
        for j = 1:N
            if j ~= i
                G_ij = get_matrix_param('G', i, j, parameter, defaultParams, active_specs);
                tau_j = get_node_param('tau', j, parameter, defaultParams, active_specs);
                Vstar_j = get_node_param('Vstar', j, parameter, defaultParams, active_specs);
                sj = 1 ./ (1 + exp(-tau_j .* (V(:, j) - Vstar_j)));
                I_cpl = I_cpl + G_ij .* sj .* (Esyn_i - V(:, i));
            end
        end
end

I_ion = gNa_i .* m(:, i).^3 .* h(:, i) .* (V(:, i) - ENa_i) ...
      + gK_i .* n(:, i).^4 .* (V(:, i) - EK_i) ...
      + gL_i .* (V(:, i) - EL_i);

value = (I_i - I_ion + I_cpl) ./ C_i;
end

function value = hh_m_rhs(variable, parameter, defaultParams, active_specs, N, i)
V = variable(:, 1:N);
[am, bm, ~, ~, ~, ~] = get_rates_scalar(V(:, i));
value = am .* (1 - variable(:, N + i)) - bm .* variable(:, N + i);
end

function value = hh_h_rhs(variable, parameter, defaultParams, active_specs, N, i)
V = variable(:, 1:N);
[~, ~, ah, bh, ~, ~] = get_rates_scalar(V(:, i));
value = ah .* (1 - variable(:, 2 * N + i)) - bh .* variable(:, 2 * N + i);
end

function value = hh_n_rhs(variable, parameter, defaultParams, active_specs, N, i)
V = variable(:, 1:N);
[~, ~, ~, ~, an, bn] = get_rates_scalar(V(:, i));
value = an .* (1 - variable(:, 3 * N + i)) - bn .* variable(:, 3 * N + i);
end

function specs = build_active_specs(active_param_names)
specs = repmat(struct('kind', '', 'base', '', 'i', [], 'j', []), 1, numel(active_param_names));
for k = 1:numel(active_param_names)
    name = char(active_param_names{k});
    parts = split(string(name), '_');
    if parts(1) == "G"
        specs(k).kind = 'matrix';
        specs(k).base = 'G';
        specs(k).i = str2double(parts(2));
        specs(k).j = str2double(parts(3));
    else
        specs(k).kind = 'node';
        specs(k).base = char(parts(1));
        specs(k).i = str2double(parts(2));
        specs(k).j = [];
    end
end
end

function params = build_active_param_vector(defaultParams, active_specs)
params = zeros(1, numel(active_specs));
for k = 1:numel(active_specs)
    spec = active_specs(k);
    if strcmp(spec.kind, 'matrix')
        params(k) = defaultParams.(spec.base)(spec.i, spec.j);
    else
        params(k) = defaultParams.(spec.base)(spec.i);
    end
end
end

function currentParams = apply_active_params(defaultParams, active_specs, parameter)
currentParams = defaultParams;
for k = 1:numel(active_specs)
    spec = active_specs(k);
    if strcmp(spec.kind, 'matrix')
        currentParams.(spec.base)(spec.i, spec.j) = parameter(k);
    else
        currentParams.(spec.base)(spec.i) = parameter(k);
    end
end
end

function value = get_node_param(base, i, parameter, defaultParams, active_specs)
value = defaultParams.(base)(i);
for k = 1:numel(active_specs)
    spec = active_specs(k);
    if strcmp(spec.kind, 'node') && strcmp(spec.base, base) && spec.i == i
        value = parameter(k);
        return
    end
end
end

function value = get_matrix_param(base, i, j, parameter, defaultParams, active_specs)
value = defaultParams.(base)(i, j);
for k = 1:numel(active_specs)
    spec = active_specs(k);
    if strcmp(spec.kind, 'matrix') && strcmp(spec.base, base) && spec.i == i && spec.j == j
        value = parameter(k);
        return
    end
end
end

function p = build_hh_p(allParams)
p = struct();
p.couplingType = allParams.couplingType;
p.C = allParams.C;
p.gNa = allParams.gNa;
p.ENa = allParams.ENa;
p.gK = allParams.gK;
p.EK = allParams.EK;
p.gL = allParams.gL;
p.EL = allParams.EL;
p.tau = allParams.tau;
p.Vstar = allParams.Vstar;
p.Esyn = allParams.Esyn;
p.G = allParams.G;
end

function values = extract_target_feature_values(items_per, features, params)
values = zeros(1, numel(items_per));
for k = 1:numel(items_per)
    item = items_per(k);
    switch item.prop
        case 'varAmp'
            values(k) = features.state.amplitude(item.idx);
        case 'varMax'
            values(k) = features.state.max(item.idx);
        case 'varMin'
            values(k) = features.state.min(item.idx);
        case 'p_Psi'
            values(k) = features.period / (2 * pi);
        case 'params'
            values(k) = params(item.idx);
    end
end
end

function name = target_feature_name(item)
name = sprintf('%s_%d', item.prop, item.idx);
end

function [am, bm, ah, bh, an, bn] = get_rates_scalar(V)
am = alpha_m_scalar(V);
bm = beta_m_scalar(V);
ah = 0.07 .* exp(-(V + 65) ./ 20);
bh = 1 ./ (1 + exp(-(V + 35) ./ 10));
an = alpha_n_scalar(V);
bn = 0.125 .* exp(-(V + 65) ./ 80);
end

function value = alpha_m_scalar(V)
value = 0.1 .* (V + 40) ./ (1 - exp(-(V + 40) ./ 10));
if isnumeric(V)
    value(isnan(value)) = 1;
end
end

function value = beta_m_scalar(V)
value = 4 .* exp(-(V + 65) ./ 18);
end

function value = alpha_n_scalar(V)
value = 0.01 .* (V + 55) ./ (1 - exp(-(V + 55) ./ 10));
if isnumeric(V)
    value(isnan(value)) = 0.1;
end
end
