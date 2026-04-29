clear
clc

%% User settings
script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, 'data');

baseline_dir = ...
    'base_N20_M10';
sweep_dir = 'truncation_order_sweep';

runtime_field = 'fmamTime';
only_plot_success = true;

data_sweep_dir = fullfile(data_dir, baseline_dir, sweep_dir);
figure_dir = fullfile(script_dir, 'figures', baseline_dir, sweep_dir);

%% Load one-result-per-file sweep data
mat_files = dir(fullfile(data_sweep_dir, '*.mat'));
if isempty(mat_files)
    error('fhnComplexity:NoPlotData', ...
        'No MAT files found in %s.', data_sweep_dir);
end

result_cells = cell(1, numel(mat_files));
for i = 1:numel(mat_files)
    data = load(fullfile(mat_files(i).folder, mat_files(i).name), 'result');
    result_cells{i} = data.result;
end
results = [result_cells{:}];

if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end

%% Draw figure for this sweep directory
[sweep_name, x_field, plot_title] = infer_sweep_plot(sweep_dir);
if strcmp(sweep_name, 'network_coupling')
    plot_network_coupling(results, figure_dir, runtime_field, only_plot_success);
else
    plot_one_sweep(results, figure_dir, sweep_name, x_field, runtime_field, ...
        plot_title, only_plot_success);
end

fprintf('Saved figures to: %s\n', figure_dir);

function [sweep_name, x_field, plot_title] = infer_sweep_plot(sweep_dir)
switch sweep_dir
    case 'sys_dim_sweep'
        sweep_name = 'dimension';
        x_field = 'N';
        plot_title = 'FMAM runtime vs system dimension';
    case 'truncation_order_sweep'
        sweep_name = 'truncation';
        x_field = 'M';
        plot_title = 'FMAM runtime vs Fourier truncation order';
    case 'target_num_sweep'
        sweep_name = 'targets';
        x_field = 'targetCount';
        plot_title = 'FMAM runtime vs target count';
    case 'err_bound_sweep'
        sweep_name = 'errBound';
        x_field = 'errBound';
        plot_title = 'FMAM runtime vs error bound';
    case 'dlambda_cap_sweep'
        sweep_name = 'dlambdaCap';
        x_field = 'dlambdaCap';
        plot_title = 'FMAM runtime vs dlambda cap';
    case 'network_coupling_sweep'
        sweep_name = 'network_coupling';
        x_field = '';
        plot_title = '';
    otherwise
        error('fhnComplexity:UnknownSweepDirectory', ...
            'Unsupported sweep directory: %s.', sweep_dir);
end
end

function plot_one_sweep(results, figure_dir, sweep_name, x_field, y_field, plot_title, only_plot_success)
mask = strcmp({results.sweepName}, sweep_name);
if only_plot_success
    mask = mask & strcmp({results.status}, 'success');
end
if ~any(mask)
    fprintf('No rows to plot for sweep: %s\n', sweep_name);
    return
end

subset = results(mask);
x = [subset.(x_field)];
y = [subset.(y_field)];
[x, order] = sort(x);
y = y(order);

fig = figure;
plot(x, y, 'o-', 'LineWidth', 1.5);
grid on
xlabel(x_field);
ylabel(y_field);
title(plot_title);
saveas(fig, fullfile(figure_dir, sprintf('fhn_%s_%s.png', sweep_name, y_field)));
end

function plot_network_coupling(results, figure_dir, y_field, only_plot_success)
mask = strcmp({results.sweepName}, 'network_coupling');
if only_plot_success
    mask = mask & strcmp({results.status}, 'success');
end
if ~any(mask)
    fprintf('No rows to plot for sweep: network_coupling\n');
    return
end

subset = results(mask);
labels = strings(1, numel(subset));
y = zeros(1, numel(subset));
for i = 1:numel(subset)
    labels(i) = sprintf('%s-%s', subset(i).netName, subset(i).couplingType);
    y(i) = subset(i).(y_field);
end

fig = figure;
bar(categorical(labels), y);
grid on
ylabel(y_field);
title('FMAM runtime by network and coupling');
saveas(fig, fullfile(figure_dir, sprintf('fhn_network_coupling_%s.png', y_field)));
end
