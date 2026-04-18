clear
clc

%% pre-setting
index_PV = 1;
period_multiplier = 1.5; % 1, 1.5, 2
unit = 'voltage'; % 'voltage' or 'concentration'

period_num = 6;
linewidth = 1;
ptNumAll = 90;
ptNum = max(8, round(ptNumAll / period_num * period_multiplier));
mksize = 25;

startidx = 50;
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
targetData = load(targetFile);
baselineData = load(baselineFile);

figure

%% ODE plot
t = targetData.TS{1};
y = targetData.TS{2} / scaleFactor;
tRel = t(:) - t(1);
Tstep = tRel(end);
tLong = tRel;
yLong = y;
for k = 1:5
    tLong = [tLong; tRel(2:end) + k * Tstep];
    yLong = [yLong; y(2:end, :)];
end
ix = find(islocalmax(yLong(:, 1)));
i0 = ix(end - 1);
i1 = ix(end);
t = tLong(i0:i1) - tLong(i0);
y = yLong(i0:i1, :);

nPoint = numel(t);
t_plot = zeros(nPoint * period_num, 1);
y_plot = zeros(nPoint * period_num, size(y, 2));
timeStart = 0;
for k = 1:period_num
    idx = (1:nPoint) + (k - 1) * nPoint;
    t_plot(idx) = t + timeStart;
    y_plot(idx, :) = y;
    timeStart = t_plot(idx(end));
end

hold on
plt_l = plot(t_plot, y_plot(:, index_PV), 'LineStyle', '-', 'LineWidth', linewidth);
grid on
box on

%% circuit plot
switch period_multiplier
    case 1
        circuitFile = fullfile(circuitDataDir, '1x_10_20.txt');
    case 1.5
        circuitFile = fullfile(circuitDataDir, '1.5x_10_20.txt');
    case 2
        circuitFile = fullfile(circuitDataDir, '2x_10_20.txt');
end
tab = readtable(circuitFile, 'VariableNamingRule', 'preserve');
tab.Properties.VariableNames{1} = 't';

t = tab.t(startidx:end) * 1e3;
V1 = tab.VF1(startidx:end) / scaleFactor;
V2 = tab.VF2(startidx:end) / scaleFactor;
V3 = tab.VF3(startidx:end) / scaleFactor;
y = [V1, V2, V3];

ix = find(islocalmax(y(:, 1)));
i0 = ix(end - 1);
i1 = ix(end);
t = t(i0:i1) - t(i0);
y = y(i0:i1, :);

tend = t(end);
tstart = t(1);
tt = (0:ptNum)' / ptNum * (tend - tstart) + tstart;
yy = zeros(numel(tt), size(y, 2));
for i = 1:size(y, 2)
    yy(:, i) = spline(t, y(:, i), tt);
end

nPoint = numel(tt);
t_plot = zeros(nPoint * period_num, 1);
y_plot = zeros(nPoint * period_num, size(yy, 2));
timeStart = 0;
for k = 1:period_num
    idx = (1:nPoint) + (k - 1) * nPoint;
    t_plot(idx) = tt + timeStart;
    y_plot(idx, :) = yy;
    timeStart = t_plot(idx(end));
end
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
