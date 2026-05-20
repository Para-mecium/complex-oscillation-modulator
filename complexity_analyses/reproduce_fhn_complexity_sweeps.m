clear
clc

script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, 'data', 'base_N20_M10');
figure_dir = fullfile(script_dir, 'figures', 'base_N20_M10', 'paper_sweeps');
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end


%% Fig. S11: system dimension
files = dir(fullfile(data_dir, 'sys_dim_sweep', '*.mat'));
x = zeros(1, numel(files));
runtime = zeros(1, numel(files));
runtime_per_iter = zeros(1, numel(files));
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name), 'result');
    x(i) = 2 * data.result.N;
    runtime(i) = data.result.fmamTime;
    runtime_per_iter(i) = data.result.fmamTime / data.result.newtonIterationsTotal;
end
[x, order] = sort(x);
runtime = runtime(order);
runtime_per_iter = runtime_per_iter(order);

fig = figure('Position', [100, 100, 620, 230]);
subplot(1, 2, 1)
h = plot(x, runtime, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'YScale', 'log')
xlim([0, 210])
ylim([1, 400])
xticks([0, 50, 100, 150, 200])
xlabel('2{\itN}', 'FontName', 'Arial')
ylabel('Runtime (s)', 'FontName', 'Arial')
subplot(1, 2, 2)
h = plot(x, runtime_per_iter, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'YScale', 'log')
xlim([0, 210])
ylim([0.04, 7])
xticks([0, 50, 100, 150, 200])
xlabel('2{\itN}', 'FontName', 'Arial')
ylabel('Avg. runtime/iter. (s)', 'FontName', 'Arial')
exportgraphics(fig, fullfile(figure_dir, 'fig_S11_system_dimension.png'), 'Resolution', 300)

%% Fig. S12: Fourier truncation order
files = dir(fullfile(data_dir, 'truncation_order_sweep', '*.mat'));
x = zeros(1, numel(files));
runtime = zeros(1, numel(files));
runtime_per_iter = zeros(1, numel(files));
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name), 'result');
    x(i) = data.result.M;
    runtime(i) = data.result.fmamTime;
    runtime_per_iter(i) = data.result.fmamTime / data.result.newtonIterationsTotal;
end
[x, order] = sort(x);
runtime = runtime(order);
runtime_per_iter = runtime_per_iter(order);

fig = figure('Position', [100, 100, 620, 230]);
subplot(1, 2, 1)
h = plot(x, runtime, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'YScale', 'log')
xlim([0, 55])
ylim([10, 80])
xticks([10, 20, 30, 40, 50])
xlabel('{\itM}', 'FontName', 'Arial')
ylabel('Runtime (s)', 'FontName', 'Arial')
subplot(1, 2, 2)
h = plot(x, runtime_per_iter, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'YScale', 'log')
xlim([0, 55])
ylim([0.1, 2])
xticks([10, 20, 30, 40, 50])
xlabel('{\itM}', 'FontName', 'Arial')
ylabel('Avg. runtime/iter. (s)', 'FontName', 'Arial')
exportgraphics(fig, fullfile(figure_dir, 'fig_S12_truncation_order.png'), 'Resolution', 300)

%% Fig. S13: number of targets
files = dir(fullfile(data_dir, 'target_num_sweep', '*.mat'));
x = zeros(1, numel(files));
runtime = zeros(1, numel(files));
runtime_per_iter = zeros(1, numel(files));
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name), 'result');
    x(i) = data.result.targetCount;
    runtime(i) = data.result.fmamTime;
    runtime_per_iter(i) = data.result.fmamTime / data.result.newtonIterationsTotal;
end
[x, order] = sort(x);
runtime = runtime(order);
runtime_per_iter = runtime_per_iter(order);

fig = figure('Position', [100, 100, 620, 230]);
subplot(1, 2, 1)
h = plot(x, runtime, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
xlim([0, 21])
ylim([10, 16])
xticks([0, 5, 10, 15, 20])
xlabel('m', 'FontName', 'Arial')
ylabel('Runtime (s)', 'FontName', 'Arial')
subplot(1, 2, 2)
h = plot(x, runtime_per_iter, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
xlim([0, 21])
ylim([0.2, 0.4])
xticks([0, 5, 10, 15, 20])
xlabel('{\itm}', 'FontName', 'Arial')
ylabel('Avg. runtime/iter. (s)', 'FontName', 'Arial')
exportgraphics(fig, fullfile(figure_dir, 'fig_S13_target_number.png'), 'Resolution', 300)

%% Fig. S14: Newton tolerance
files = dir(fullfile(data_dir, 'err_bound_sweep', '*.mat'));
x = zeros(1, numel(files));
runtime = zeros(1, numel(files));
runtime_per_iter = zeros(1, numel(files));
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name), 'result');
    x(i) = data.result.errBound;
    runtime(i) = data.result.fmamTime;
    runtime_per_iter(i) = data.result.fmamTime / data.result.newtonIterationsTotal;
end
[x, order] = sort(x);
runtime = runtime(order);
runtime_per_iter = runtime_per_iter(order);

fig = figure('Position', [100, 100, 620, 230]);
subplot(1, 2, 1)
h = plot(x, runtime, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'XScale', 'log')
xlim([1e-10, 1e-4])
ylim([10, 16])
xlabel('\epsilon', 'FontName', 'Arial')
ylabel('Runtime (s)', 'FontName', 'Arial')
subplot(1, 2, 2)
h = plot(x, runtime_per_iter, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'XScale', 'log')
xlim([1e-10, 1e-4])
ylim([0.2, 0.4])
xlabel('\epsilon', 'FontName', 'Arial')
ylabel('Avg. runtime/iter. (s)', 'FontName', 'Arial')
exportgraphics(fig, fullfile(figure_dir, 'fig_S14_newton_tolerance.png'), 'Resolution', 300)

%% Fig. S15: continuation step cap
files = dir(fullfile(data_dir, 'dlambda_cap_sweep', '*.mat'));
x = zeros(1, numel(files));
runtime = zeros(1, numel(files));
runtime_per_iter = zeros(1, numel(files));
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name), 'result');
    x(i) = round(1 / data.result.dlambdaCap);
    runtime(i) = data.result.fmamTime;
    runtime_per_iter(i) = data.result.fmamTime / data.result.newtonIterationsTotal;
end
[x, order] = sort(x);
runtime = runtime(order);
runtime_per_iter = runtime_per_iter(order);

fig = figure('Position', [100, 100, 620, 230]);
subplot(1, 2, 1)
h = plot(x, runtime, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
set(gca, 'YScale', 'log')
xlim([0, 210])
ylim([6, 200])
xticks([0, 50, 100, 150, 200])
xlabel('{\itK}', 'FontName', 'Arial')
ylabel('Runtime (s)', 'FontName', 'Arial')
subplot(1, 2, 2)
h = plot(x, runtime_per_iter, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6);
h.MarkerFaceColor = h.Color;
h.MarkerEdgeColor = 'w';
grid on
box on
xlim([0, 210])
ylim([0.2, 0.4])
xticks([0, 50, 100, 150, 200])
xlabel('{\itK}', 'FontName', 'Arial')
ylabel('Avg. runtime/iter. (s)', 'FontName', 'Arial')
exportgraphics(fig, fullfile(figure_dir, 'fig_S15_step_cap.png'), 'Resolution', 300)
