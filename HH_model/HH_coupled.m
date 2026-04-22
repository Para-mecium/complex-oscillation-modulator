% Coupled Hodgkin-Huxley model:
% compare slow external-current change and sudden external-current drop.

clear; clc; 

%% 1. Local configuration
N = 2;

% random seed
seed = 1;
rng(seed, 'twister');

% Shared HH parameters written in vector form.
C0 = 1.0;
gNa0 = 120.0; ENa0 = 50.0;
gK0 = 36.0;   EK0 = -77.0;
gL0 = 0.3;    EL0 = -54.387;

% Current levels (identical across neurons, kept in vector form).
I_high_scalar = 120;
I_low_scalar = 80;

% Shared initial condition, replicated to all neurons.
V0 = 0.2;
m0 = 0.5;
h0 = 0.2;
n0 = 0.5;

% Solver options.
RelTol = 1e-6;
AbsTol = 1e-6;

% Time settings
T_start = 10000;
T_slow_total = 2*T_start;
T_step_total = 2*T_start;
t_hold = T_start;
t_ramp_end = 1.5*T_start;

t_init_end = T_start;
t_search_start = T_start;
t_search_end = T_start + 15;
n_search = 300;
t_test_duration = 50;

% Observation settings.
obs = @(Y) mean(Y(:, 1:N), 2);
obsLabel = 'Mean Voltage (mV)';
V_tol = -60;

% Viz settings
t_TS_plot_span = 2*T_start;


%% 2. Build parameters
p = struct();

G_scale = 0.01;
if N > 1
    p.G = G_scale * (2*rand(N,N)-1);
end

p.C = C0 * ones(N, 1);
p.gNa = gNa0 * ones(N, 1);
p.ENa = ENa0 * ones(N, 1);
p.gK = gK0 * ones(N, 1);
p.EK = EK0 * ones(N, 1);
p.gL = gL0 * ones(N, 1);
p.EL = EL0 * ones(N, 1);

I_high = [120;120];
I_low = [9.40124960614535;	5.99661111564136]; %

y0 = [ ...
    V0 * randn(N, 1); ...
    m0 * ones(N, 1); ...
    h0 * ones(N, 1); ...
    n0 * ones(N, 1)];

options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

%% 4. Slow current decrease (control)
I_slow = @(t) linear_current_profile(t, I_high, I_low, t_hold, t_ramp_end);

fprintf('正在计算缓慢变化情形...\n');
[t_slow, y_slow] = ode15s(@(t, y) build_model(t, y, I_slow(t), p, N), ...
    [0, T_slow_total], y0, options);

obs_slow = obs(y_slow);
% obs_slow = y_slow(:,1) - y_slow(:,2);
I_slow_val = mean(I_slow(t_slow), 2);

% %% 5. Sudden current drop: search for a suitable switching time
% fprintf('正在搜索突降时刻...\n');
% [t_init, y_init] = ode15s(@(t, y) build_model(t, y, I_high, p, N), ...
%     [0, t_init_end], y0, options);
% y_steady = y_init(end, :).';
% 
% t_search = linspace(t_search_start, t_search_end, n_search);
% [~, y_search] = ode15s(@(t, y) build_model(t, y, I_high, p, N), ...
%     t_search, y_steady, options);
% 
% t_drop = t_search_start;
% found = false;
% for i = 1:numel(t_search)
%     y_test_start = y_search(i, :).';
%     [~, y_test] = ode15s(@(t, y) build_model(t, y, I_low, p, N), ...
%         [0, t_test_duration], y_test_start, options);
% 
%     obs_test = obs(y_test);
%     tail_idx = floor(numel(obs_test) / 2) + 1:numel(obs_test);
%     obs_tail = obs_test(tail_idx);
% 
%     if max(obs_tail) < V_tol
%         t_drop = t_search(i);
%         found = true;
%         break;
%     end
% end

% if found
%     fprintf('找到合适的突降时刻: t_drop = %.3f ms\n', t_drop);
% else
%     fprintf('警告: 未找到满足振幅阈值的突降时刻，使用默认值。\n');
% end

t_drop = t_search_start;
I_step = @(t) step_current_profile(t, I_high, I_low, t_drop);

[t_step, y_step] = ode15s(@(t, y) build_model(t, y, I_step(t), p, N), ...
    [0, T_step_total], y0, options);

obs_step = obs(y_step);
I_step_val = mean(I_step(t_step), 2);

%% 6. Visualization
figure('Position', [50, 100, 1200, 600], ...
    'Name', 'Coupled Hodgkin-Huxley: Slow vs Sudden Current Change');

I_ylim = [min(0, min(I_low) - 10), max(I_high) + 10];

subplot(2, 2, 1);
plot(t_slow, obs_slow, 'b-', 'LineWidth', 1.2);
title('Slow Current Decrease', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(obsLabel, 'FontSize', 11);
grid on;
xlim([0, T_slow_total]);

subplot(2, 2, 3);
plot(t_slow, I_slow_val, 'r-', 'LineWidth', 2);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Mean Current (\muA/cm^2)', 'FontSize', 11);
ylim(I_ylim);
grid on;
xlim([0, T_slow_total]);

subplot(2, 2, 2);
plot(t_step, obs_step, 'b-', 'LineWidth', 1.2);
title('Sudden Current Drop', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(obsLabel, 'FontSize', 11);
grid on;
xlim([0, T_step_total]);
xline(t_drop, 'k--', sprintf('Drop to I=%g', I_low_scalar), ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 10);

subplot(2, 2, 4);
plot(t_step, I_step_val, 'r-', 'LineWidth', 2);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Mean Current (\muA/cm^2)', 'FontSize', 11);
ylim(I_ylim);
grid on;
xlim([0, T_step_total]);

sgtitle(sprintf('Coupled HH Model (N = %d)', N), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% 7. Voltage traces of all neurons
V_slow = y_slow(:, 1:N);
V_step = y_step(:, 1:N);
colors = lines(N);

slow_plot_idx = t_slow > (t_slow-t_TS_plot_span);
step_plot_idx = t_step > (t_step-t_TS_plot_span);
t_slow_TS_plot_start = T_slow_total - t_TS_plot_span;
t_step_TS_plot_start = T_step_total - t_TS_plot_span;

figure('Position', [80, 120, 1200, 500], ...
    'Name', 'Neuron Voltage Traces');

subplot(1, 2, 1);
hold on;
for i = 1:N
    plot(t_slow(slow_plot_idx)-t_slow_TS_plot_start, V_slow(slow_plot_idx, i), 'LineWidth', 1.2, 'Color', colors(i, :));
end
hold off;
title('Slow Current Decrease: All Voltages', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Voltage (mV)', 'FontSize', 11);
grid on;
xlim([0, t_TS_plot_span]);

subplot(1, 2, 2);
hold on;
for i = 1:N
    plot(t_step(step_plot_idx)-t_step_TS_plot_start, V_step(step_plot_idx, i), 'LineWidth', 1.2, 'Color', colors(i, :));
end
hold off;
title('Sudden Current Drop: All Voltages', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Voltage (mV)', 'FontSize', 11);
grid on;
xlim([0, t_TS_plot_span]);
xline(t_drop, 'k--', 'LineWidth', 1);

%% ================= Local helpers ================= %%

function I_t = linear_current_profile(t, I_high, I_low, t_hold, t_ramp_end)
t = t(:);
n_times = numel(t);
N = numel(I_high);
I_t = zeros(n_times, N);

for k = 1:n_times
    tk = t(k);
    if tk < t_hold
        I_t(k, :) = I_high.';
    elseif tk <= t_ramp_end
        lambda = (tk - t_hold) / (t_ramp_end - t_hold);
        I_t(k, :) = ((1 - lambda) * I_high + lambda * I_low).';
    else
        I_t(k, :) = I_low.';
    end
end

if n_times == 1
    I_t = I_t.';
end
end

function I_t = step_current_profile(t, I_high, I_low, t_drop)
t = t(:);
n_times = numel(t);
N = numel(I_high);
I_t = zeros(n_times, N);

for k = 1:n_times
    if t(k) < t_drop
        I_t(k, :) = I_high.';
    else
        I_t(k, :) = I_low.';
    end
end

if n_times == 1
    I_t = I_t.';
end
end
