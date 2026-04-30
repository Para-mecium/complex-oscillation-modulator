clear
clc

dynamicName = 'GRN';
netName = 'ER';
N = 100;
n_per = 50;
enableCustomPlotSettings = false;
color_ori = [0, 0.4470, 0.7410];
color_per = [0.8500, 0.3250, 0.0980];

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
if ~isempty(plotSettings.axisRange)
    axis(plotSettings.axisRange)
end
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
if ~isempty(plotSettings.axisRange)
    axis(plotSettings.axisRange)
end
xlabel(plotSettings.xname)
ylabel(plotSettings.yname)
set(gca, 'fontsize', 8)
format_figure(gcf)
print(gcf, fullfile(figDir, ...
    ['Phase_' dynamicName '_n_per = ' num2str(n_per) '_' netName '.pdf']), '-dpdf');

figure
hold on
for i = 1:N
    patch([i - 0.4, i + 0.4, i + 0.4, i - 0.4], [0, 0, amp_v(i), amp_v(i)], ...
        color_per, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
for i = 1:N
    patch([i - 0.4, i + 0.4, i + 0.4, i - 0.4], [0, 0, amp_v_origin(i), amp_v_origin(i)], ...
        color_ori, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
xlim([0, N + 1]);
box on
xlabel('Node index')
ylabel(plotSettings.amplitudeLabel)
set(gca, 'fontsize', 8)
format_figure(gcf)
print(gcf, fullfile(figDir, ...
    ['Bar_' dynamicName '_n_per = ' num2str(n_per) '_' netName '.pdf']), '-dpdf');

figure
hold on
if isempty(plotSettings.histogramEdges)
    histogram(amp_v, FaceColor=color_per, FaceAlpha=0.5, EdgeAlpha=0.5)
    histogram(amp_v_origin, FaceColor=color_ori, FaceAlpha=0.5, EdgeAlpha=0.5)
else
    histogram(amp_v, FaceColor=color_per, BinEdges=plotSettings.histogramEdges, FaceAlpha=0.5, EdgeAlpha=0.5)
    histogram(amp_v_origin, FaceColor=color_ori, BinEdges=plotSettings.histogramEdges, FaceAlpha=0.5, EdgeAlpha=0.5)
end
box on
grid on
xlabel(plotSettings.amplitudeLabel)
ylabel('Number of nodes')
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
            'histogramEdges', linspace(15, 30, 50), ...
            'histogramAxis', []);
    otherwise
        plotSettings = struct( ...
            'axisRange', [-0.2, 1.4, 0.4, 0.8], ...
            'xname', '{\it V} (a.u.)', ...
            'yname', '{\it W} (a.u.)', ...
            'amplitudeLabel', 'Amplitude of {\it V} (a.u.)', ...
            'histogramEdges', linspace(0.56, 0.8, 35), ...
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
            'histogramAxis', []);
    otherwise
        plotSettings = struct( ...
            'axisRange', [], ...
            'xname', '{\it V} (a.u.)', ...
            'yname', '{\it W} (a.u.)', ...
            'amplitudeLabel', 'Amplitude of {\it V} (a.u.)', ...
            'histogramEdges', [], ...
            'histogramAxis', []);
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
