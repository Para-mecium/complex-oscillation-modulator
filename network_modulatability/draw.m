clear
clc

dynamicName = 'GRN';
netName = 'SW';
N = 100;
n_per = 50;
enableCustomPlotSettings = false;
color_ori = [0, 0.4470, 0.7410];
color_per = [0.8500, 0.3250, 0.0980];
phasePaddingFraction = 0.05;

scriptDir = fileparts(mfilename('fullpath'));
folderName = fullfile(scriptDir, [dynamicName '_' netName]);
figDir = fullfile(scriptDir, 'temp_fig');
ensure_directory(figDir);

baseline = networkexp.load_single_result( ...
    fullfile(folderName, sprintf('TS_init_%s (%s, N = %d).mat', ...
    dynamicName, netName, N)));
perturbed = networkexp.load_single_result( ...
    fullfile(folderName, sprintf('TS_per_%d_%s (%s, N = %d).mat', ...
    n_per, dynamicName, netName, N)));

amp_v_origin = baseline.amplitude(1:N);
amp_v = perturbed.amplitude(1:N);
plotSettings = get_plot_settings(dynamicName, enableCustomPlotSettings);

t_ori = baseline.TS{1};
TSvar_ori = baseline.TS{2};
TSvar = perturbed.TS{2};
if isempty(t_ori) || isempty(TSvar_ori) || isempty(perturbed.TS{1}) || isempty(TSvar)
    error('draw:MissingPeriodicOrbit', ...
        'One of the loaded runs does not contain a periodic orbit trajectory.');
end
phaseAxisRange = get_phase_axis_range(TSvar, N, phasePaddingFraction);

fig = figure;
hold on
for i = 1:N
    V_i = TSvar_ori(:, i);
    W_i = TSvar_ori(:, i + N);
    k = boundary(V_i, W_i, 0);
    pgon = polyshape(V_i(k), W_i(k));
    plot(pgon, FaceAlpha=0, EdgeColor=color_ori, EdgeAlpha=1)
end
grid on
box on
axis(phaseAxisRange)
apply_phase_tick_settings(dynamicName)
xlabel(plotSettings.xname)
ylabel(plotSettings.yname)
set(gca, 'fontsize', 8)
format_figure(gcf)
print(gcf, fullfile(figDir, ['Phase_' dynamicName '_init_' netName '.pdf']), '-dpdf');

fig = figure;
hold on
for i = 1:N
    V_i = TSvar(:, i);
    W_i = TSvar(:, i + N);
    k = boundary(V_i, W_i, 0);
    pgon = polyshape(V_i(k), W_i(k));
    if i <= n_per
        plot(pgon, FaceAlpha=0, EdgeColor=color_per, EdgeAlpha=1)
    else
        plot(pgon, FaceAlpha=0, EdgeColor=color_ori, EdgeAlpha=1)
    end
end
grid on
box on
axis(phaseAxisRange)
apply_phase_tick_settings(dynamicName)
xlabel(plotSettings.xname)
ylabel(plotSettings.yname)
set(gca, 'fontsize', 8)
format_figure(gcf)
print(gcf, fullfile(figDir, ...
    ['Phase_' dynamicName '_n_per = ' num2str(n_per) '_' netName '.pdf']), '-dpdf');

figure
hold on
for i = 1:N
    hPerBar = patch([i - 0.4, i + 0.4, i + 0.4, i - 0.4], [0, 0, amp_v(i), amp_v(i)], ...
        color_per, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
for i = 1:N
    hOriBar = patch([i - 0.4, i + 0.4, i + 0.4, i - 0.4], [0, 0, amp_v_origin(i), amp_v_origin(i)], ...
        color_ori, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
xlim([0, N + 1]);
box on
xlabel('Node index')
ylabel(plotSettings.amplitudeLabel)
legend([hOriBar, hPerBar], {'Original', 'Modulated'}, 'Location', 'best')
set(gca, 'fontsize', 8)
format_figure(gcf)
print(gcf, fullfile(figDir, ...
    ['Bar_' dynamicName '_n_per = ' num2str(n_per) '_' netName '.pdf']), '-dpdf');

figure
hold on
histogramEdges = get_histogram_edges(amp_v, amp_v_origin, dynamicName, plotSettings);
hPerHist = histogram(amp_v, FaceColor=color_per, BinEdges=histogramEdges, FaceAlpha=0.5, EdgeAlpha=0.5);
hOriHist = histogram(amp_v_origin, FaceColor=color_ori, BinEdges=histogramEdges, FaceAlpha=0.5, EdgeAlpha=0.5);
box on
grid on
xlabel(plotSettings.amplitudeLabel)
ylabel('Number of nodes')
legend([hOriHist, hPerHist], {'Original', 'Modulated'}, 'Location', 'best')
if ~isempty(plotSettings.histogramAxis)
    axis(plotSettings.histogramAxis)
end
set(gca, 'fontsize', 8)
format_figure(gcf)
print(gcf, fullfile(figDir, ...
    ['Histo_' dynamicName '_n_per = ' num2str(n_per) '_' netName '.pdf']), '-dpdf');

function plotSettings = get_plot_settings(dynamicName, enableCustomPlotSettings)
plotSettings = get_default_plot_settings(dynamicName);
if ~enableCustomPlotSettings
    return;
end

switch upper(dynamicName)
    case 'GRN'
        plotSettings = struct( ...
            'axisRange', [-1, 55, 0, 0.5], ...
            'xname', 'Protein (a.u.)', ...
            'yname', 'mRNA (a.u.)', ...
            'amplitudeLabel', 'Amplitude of protein (a.u.)', ...
            'histogramEdges', linspace(15, 30, 61), ...
            'histogramRange', [], ...
            'histogramBinCount', 60, ...
            'histogramAxis', []);
    otherwise
        plotSettings = struct( ...
            'axisRange', [-0.2, 1.4, 0.4, 0.8], ...
            'xname', '{\it V} (a.u.)', ...
            'yname', '{\it W} (a.u.)', ...
            'amplitudeLabel', 'Amplitude of {\it V} (a.u.)', ...
            'histogramEdges', linspace(0.56, 0.8, 61), ...
            'histogramRange', [], ...
            'histogramBinCount', 60, ...
            'histogramAxis', []);
end
end

function plotSettings = get_default_plot_settings(dynamicName)
switch upper(dynamicName)
    case 'GRN'
        plotSettings = struct( ...
            'axisRange', [], ...
            'xname', 'Protein (a.u.)', ...
            'yname', 'mRNA (a.u.)', ...
            'amplitudeLabel', 'Amplitude of protein (a.u.)', ...
            'histogramEdges', [], ...
            'histogramRange', [], ...
            'histogramBinCount', 60, ...
            'histogramAxis', []);
    otherwise
        plotSettings = struct( ...
            'axisRange', [], ...
            'xname', '{\it V} (a.u.)', ...
            'yname', '{\it W} (a.u.)', ...
            'amplitudeLabel', 'Amplitude of {\it V} (a.u.)', ...
            'histogramEdges', [], ...
            'histogramRange', [], ...
            'histogramBinCount', 60, ...
            'histogramAxis', []);
end
end

function apply_phase_tick_settings(dynamicName)
if strcmpi(dynamicName, 'FHN')
    yTickValues = 0.3:0.1:0.7;
    yLimits = ylim;
    ylim([min(yLimits(1), yTickValues(1)), max(yLimits(2), yTickValues(end))])
    yticks(yTickValues)
    yticklabels(compose('%.1f', yTickValues))
end
end

function axisRange = get_phase_axis_range(TSvar, nodeCount, paddingFraction)
xValues = TSvar(:, 1:nodeCount);
yValues = TSvar(:, nodeCount + 1:2 * nodeCount);
xValues = xValues(isfinite(xValues));
yValues = yValues(isfinite(yValues));
if isempty(xValues) || isempty(yValues)
    error('draw:MissingPhaseAxisData', ...
        'Cannot build phase axis range because modulated trajectory values are empty or non-finite.');
end

xLimits = pad_limits(min(xValues), max(xValues), paddingFraction);
yLimits = pad_limits(min(yValues), max(yValues), paddingFraction);
axisRange = [xLimits, yLimits];
end

function limits = pad_limits(dataMin, dataMax, paddingFraction)
dataRange = dataMax - dataMin;
if dataRange == 0
    padding = max(abs(dataMin) * paddingFraction, eps(max(abs([dataMin, 1]))));
else
    padding = dataRange * paddingFraction;
end
limits = [dataMin - padding, dataMax + padding];
end

function histogramEdges = get_histogram_edges(amp_v, amp_v_origin, dynamicName, plotSettings)
if ~isempty(plotSettings.histogramEdges)
    histogramEdges = plotSettings.histogramEdges;
    return
end

values = [amp_v(:); amp_v_origin(:)];
values = values(isfinite(values));
if isempty(values)
    error('draw:MissingHistogramData', ...
        'Cannot build histogram edges because amplitude values are empty or non-finite.');
end

binCount = 50;
if isfield(plotSettings, 'histogramBinCount') && ~isempty(plotSettings.histogramBinCount)
    binCount = max(50, plotSettings.histogramBinCount);
end

histogramRange = [];
if isfield(plotSettings, 'histogramRange')
    histogramRange = plotSettings.histogramRange;
end
if isempty(histogramRange)
    histogramRange = get_default_histogram_range(dynamicName, values);
end

dataMin = min(values);
dataMax = max(values);
if dataMin < histogramRange(1) || dataMax > histogramRange(2)
    padding = max(0.05 * (dataMax - dataMin), eps(max(abs([dataMin, dataMax, 1]))));
    histogramRange = [min(histogramRange(1), dataMin - padding), ...
        max(histogramRange(2), dataMax + padding)];
end
if histogramRange(1) == histogramRange(2)
    padding = max(0.05 * abs(histogramRange(1)), eps(max(abs([histogramRange(1), 1]))));
    histogramRange = histogramRange + [-padding, padding];
end

histogramEdges = linspace(histogramRange(1), histogramRange(2), binCount + 1);
end

function histogramRange = get_default_histogram_range(dynamicName, values)
switch upper(dynamicName)
    case 'FHN'
        histogramRange = [0.55, 0.65];
    case 'GRN'
        if max(values) > 1
            histogramRange = [15, 30];
        else
            histogramRange = [0.20, 0.35];
        end
    otherwise
        dataMin = min(values);
        dataMax = max(values);
        padding = max(0.05 * (dataMax - dataMin), eps(max(abs([dataMin, dataMax, 1]))));
        histogramRange = [dataMin - padding, dataMax + padding];
end
end

function ensure_directory(folderPath)
if exist(folderPath, 'dir') ~= 7
    mkdir(folderPath);
end
end

function format_figure(fig)
fig.PaperUnits = "centimeters";
width = 6;
height = 6 * 0.75;
fig.PaperSize = [width, height];
fig.PaperPosition = [0, 0, width, height];
end
