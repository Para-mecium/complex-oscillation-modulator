index_PV = 1;
period_multiplier = 2; % 1, 1.5, 2
unit = 'concentration'; % 'voltage' or 'concentration'

period_num = 6;
linewidth = 1;
ptNumAll = 90;
ptNum = max(8, round(ptNumAll / period_num * period_multiplier));
mksize = 25;
startidx = 1;
scaleToAU = 0.06;

colors_l = {'#CBBBC1', '#E4B7BC', '#F5E4C8'};
colors_s = {'#551F33', '#BD4146', '#ECC68C'};

scriptDir = fileparts(mfilename('fullpath'));
circuitDataDir = fullfile(scriptDir, 'circuit_data');

addpath(scriptDir);
switch lower(unit)
    case 'voltage'
        scaleFactor = 1;
        yAxisLabel = 'Voltage (V)';
    case 'concentration'
        scaleFactor = scaleToAU;
        yAxisLabel = 'Concentration (a.u.)';
end

baselineFile = fullfile(scriptDir, 'period_target_1p0x.mat');
switch period_multiplier
    case 1
        targetFile = fullfile(scriptDir, 'period_target_1p0x.mat');
    case 1.5
        targetFile = fullfile(scriptDir, 'period_target_1p5x.mat');
    case 2
        targetFile = fullfile(scriptDir, 'period_target_2p0x.mat');
end
if ~isfile(baselineFile) || ~isfile(targetFile)
    disp('Target cache missing. Generate with Modulation_bi_target...');
    Modulation_bi_target();
end

targetData = load(targetFile);
baselineData = load(baselineFile);

figure

%% ODE plot
t = targetData.TS{1};
y = targetData.TS{2} / scaleFactor;
[t_plot, y_plot] = repeat_periodic_trajectory(t, y, period_num);

hold on
plt_l = plot(t_plot, y_plot(:, index_PV), 'LineStyle', '-', 'LineWidth', linewidth);
grid on
box on

%% circuit plot
switch period_multiplier
    case 1
        circuitFile = fullfile(circuitDataDir, '1x_10.txt');
    case 1.5
        circuitFile = fullfile(circuitDataDir, '1.5x_10.txt');
    case 2
        circuitFile = fullfile(circuitDataDir, '2x_10.txt');
end
tab = readtable(circuitFile, 'VariableNamingRule', 'preserve');
tab.Properties.VariableNames{1} = 't';

t = tab.t(startidx:end) * 1e3;
V1 = tab.VF1(startidx:end) / scaleFactor;
V2 = tab.VF2(startidx:end) / scaleFactor;
V3 = tab.VF3(startidx:end) / scaleFactor;
y = [V1, V2, V3];

localmax = islocalmax(y(:, 1));
index_max = find(localmax);
if numel(index_max) < 2
    error('draw_TS:InsufficientCircuitPeaks', ...
        'At least two local maxima of V1 are required to extract one period.');
end

index_start = index_max(end - 1);
index_end = index_max(end);

t = t(index_start:index_end) - t(index_start);
y = y(index_start:index_end, :);

tend = t(end);
tstart = t(1);
tt = (0:ptNum)' / ptNum * (tend - tstart) + tstart;
yy = zeros(numel(tt), size(y, 2));
for i = 1:size(y, 2)
    yy(:, i) = spline(t, y(:, i), tt);
end

[t_plot, y_plot] = repeat_periodic_trajectory(tt, yy, period_num);
plt_s = scatter(t_plot, y_plot(:, index_PV), mksize, 'Marker', 'x');

%% figure properties
T_0 = baselineData.period;

xlabel('Time (a.u.)', 'FontName', 'Arial')
ylabel(yAxisLabel, 'FontName', 'Arial')

XTick_2DTS = 0:T_0:(period_num * T_0);
XTickLabel_2DTS = repmat({''}, 1, numel(XTick_2DTS));
XTickLabel_2DTS{1} = '0';
for tickIdx = 2:numel(XTick_2DTS)
    periodCount = tickIdx - 1;
    if mod(periodCount, 2) == 0
        XTickLabel_2DTS{tickIdx} = [num2str(periodCount), 'T_0'];
    end
end

switch index_PV
    case 1
        YTick_2DTS = [0, 5, 10, 15, 20, 25];
    case {2, 3}
        YTick_2DTS = [0, 15, 30, 45, 60, 75];
end
switch lower(unit)
    case 'voltage'
        YTick_2DTS = YTick_2DTS * scaleToAU;
end
YTickLabel_2DTS = cellstr(compose('%.3g', YTick_2DTS));
axis([0, T_0 * period_num, 0, YTick_2DTS(end)])

set(gca, 'FontSize', 10, 'XTick', XTick_2DTS, 'YTick', YTick_2DTS, ...
    'XTickLabel', XTickLabel_2DTS, 'YTickLabel', YTickLabel_2DTS)
set(plt_l, 'Color', colors_l{index_PV})
set(plt_s, 'MarkerFaceColor', colors_s{index_PV}, 'MarkerEdgeColor', colors_s{index_PV})

%% local functions
function [tPlot, yPlot] = repeat_periodic_trajectory(t, y, repeatCount)
t = t(:);
y = double(y);
nPoint = numel(t);
nVar = size(y, 2);

tPlot = zeros(nPoint * repeatCount, 1);
yPlot = zeros(nPoint * repeatCount, nVar);
timeStart = 0;

for k = 1:repeatCount
    idx = (1:nPoint) + (k - 1) * nPoint;
    tPlot(idx) = t + timeStart;
    yPlot(idx, :) = y;
    timeStart = tPlot(idx(end));
end
end
