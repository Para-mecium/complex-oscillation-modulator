clear
clc

cls = 'learned'; % 'learned' or 'init'
unit = 'voltage'; % 'voltage' or 'concentration'

%% pre-setting
lineWidth = 1;
index_PV = 1;
mksize = 25;
ptNum = 200;
stepsize = 1;
scaleToAU = 0.06;

colors_l = {'#CBBBC1','#E4B7BC','#F5E4C8'};
colors_s = {'#551F33','#BD4146','#ECC68C'};

scriptDir = fileparts(mfilename('fullpath'));
circuitDir = fileparts(scriptDir);

switch lower(unit)
    case 'voltage'
        scaleFactor = 1;
        yAxisLabel = 'Voltage (V)';
    case 'concentration'
        scaleFactor = scaleToAU;
        yAxisLabel = 'Concentration (a.u.)';
end

%% load ODE trajectory
switch lower(cls)
    case 'learned'
        learnedData = load(fullfile(circuitDir, "learnedData_ODE.mat"));
        t_ode = learnedData.TS{1};
        y_ode = learnedData.TS{2};
        repeatCount = 2;

    case 'init'
        initDataODE = load(fullfile(scriptDir, "initData_ODE.mat"));
        t_ode = initDataODE.TS{1};
        y_ode = initDataODE.TS{2};
        repeatCount = 6;
end

%% load circuit reference trajectory
initDataCircuit = load(fullfile(scriptDir, "initData_circuit.mat"));
t_circuit = initDataCircuit.t;
y_circuit = initDataCircuit.y;

%% ODE curve
[t_plot_ode, y_plot_ode] = build_periodic_curve(t_ode, y_ode, index_PV, repeatCount, scaleFactor);

figure
hold on
plt_l = plot(t_plot_ode, y_plot_ode, 'LineStyle', '-', 'LineWidth', lineWidth);
grid on
box on

%% circuit scatter
[tDense, yDense] = build_periodic_curve(t_circuit, y_circuit, index_PV, 2, scaleFactor);
tEnd = tDense(end);
t_plot_circuit = (0:ptNum-1)' / ptNum * tEnd;
y_plot_circuit = zeros(ptNum, size(yDense, 2));
for i = 1:size(yDense, 2)
    y_plot_circuit(:, i) = spline(tDense, yDense(:, i), t_plot_circuit);
end

plt_s = scatter(t_plot_circuit(1:stepsize:end), y_plot_circuit(1:stepsize:end, :), ...
    mksize, 'Marker', 'x');

xlabel('Time (a.u.)', 'FontName', 'Arial')
ylabel(yAxisLabel, 'FontName', 'Arial')

for i = 1:3
    set(plt_l(i), 'Color', colors_l{i});
    set(plt_s(i), 'MarkerFaceColor', colors_s{i}, 'MarkerEdgeColor', colors_s{i});
end

xMax = max([t_plot_ode(:); t_plot_circuit(:)]);
yMin = min([y_plot_ode(:); y_plot_circuit(:)]);
yMax = max([y_plot_ode(:); y_plot_circuit(:)]);
yPad = max(1e-6, 0.05 * (yMax - yMin));

yLower = max(0, yMin - yPad);
yUpper = yMax + yPad;
axis([0, xMax, yLower, yUpper]);

xTicks = linspace(0, xMax, 5);
yTicks = linspace(yLower, yUpper, 6);
set(gca, 'FontSize', 10, 'XTick', xTicks, 'YTick', yTicks, ...
    'XTickLabel', compose('%.2g', xTicks), ...
    'YTickLabel', compose('%.2g', yTicks));

%% local functions
function [tPlot, yPlot] = build_periodic_curve(t, y, indexPV, repeatCount, scaleFactor)
    t = t(:);
    [~, idx_max] = max(y(:, indexPV));

    y = [y(idx_max:end-1, :); y(1:idx_max, :)];
    y = y / scaleFactor;

    t = [t(idx_max:end-1) - t(idx_max); t(1:idx_max) + t(end) - t(idx_max)];

    tPlot = [];
    yPlot = [];
    tStart = 0;
    for k = 1:repeatCount
        if k == 1
            tSeg = t + tStart;
            ySeg = y;
        else
            tSeg = t(2:end) + tStart;
            ySeg = y(2:end, :);
        end
        tPlot = [tPlot; tSeg]; %#ok<AGROW>
        yPlot = [yPlot; ySeg]; %#ok<AGROW>
        tStart = tPlot(end);
    end
end
