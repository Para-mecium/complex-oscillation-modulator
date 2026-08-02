clear
clc
%%
% 'parameter' or 'orbit_l2'
errorMetric = 'orbit_l2';
foldYAxis = should_fold_y_axis(errorMetric);
foldAxisPadding = 0.18;

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
resultsRootDir = fullfile(scriptDir, 'results');
scaleLevels = [0.5 0.75 1 1.25 1.5];

metricSpec = build_metric_spec(errorMetric, circuitDir);

scaleValues = zeros(0, 1);
distancesByScale = {};
statsByScale = struct('numTotal', {}, 'numSuccess', {}, ...
    'successRate', {}, 'numMetricValues', {});

for scaleIdx = 1:numel(scaleLevels)
    scale = scaleLevels(scaleIdx);
    resultScaleDir = fullfile(resultsRootDir, scale_dir_name(scale));
    statsData = load(fullfile(resultScaleDir, metricSpec.statsFileName));
    values = collect_metric_values(statsData, metricSpec);

    scaleValues(end + 1, 1) = statsData.scale; %#ok<SAGROW>
    distancesByScale{end + 1, 1} = values; %#ok<SAGROW>
    statsByScale(end + 1).numTotal = statsData.numTotal; %#ok<SAGROW>
    statsByScale(end).numSuccess = statsData.numSuccess;
    statsByScale(end).successRate = statsData.successRate;
    statsByScale(end).numMetricValues = numel(values);
end

allDistances = vertcat(distancesByScale{:});
displayFloor = choose_display_floor(allDistances, metricSpec.useLogScale);

figWidth = max(640, 150 * numel(distancesByScale));
fig = figure('Color', 'w', 'Position', [100, 100, figWidth, 520]);
if foldYAxis
    [topAx, bottomAx] = create_folded_axes(fig);
    [lowerYLim, upperYLim, lowerCut, upperCut] = folded_y_limits( ...
        allDistances, foldAxisPadding, metricSpec.useLogScale);
else
    ax = axes(fig);
    hold(ax, 'on')
    grid(ax, 'on')
    box(ax, 'on')
end

colors = lines(numel(distancesByScale));
scale15Group = scaleValues == 1.5;
if any(scale15Group)
    colors(scale15Group, :) = repmat([0.2000, 0.2000, 0.2000], nnz(scale15Group), 1);
end
rng(1, 'twister');

for scaleIdx = 1:numel(distancesByScale)
    values = distancesByScale{scaleIdx};
    plotValues = max(values, displayFloor);
    numFloorClipped = nnz(values < displayFloor);
    if foldYAxis
        draw_raincloud_group(bottomAx, scaleIdx, plotValues(plotValues <= lowerCut), ...
            colors(scaleIdx, :), metricSpec.useLogDensity);
        draw_scatter_group(topAx, scaleIdx, plotValues(plotValues >= upperCut), colors(scaleIdx, :));
    else
        draw_raincloud_group(ax, scaleIdx, plotValues, colors(scaleIdx, :), metricSpec.useLogDensity);
    end

    fprintf(['scale=%.6g: success=%d/%d (%.3f), n=%d, ' ...
        'median%s=%.6g, floorClipped=%d\n'], ...
        scaleValues(scaleIdx), ...
        statsByScale(scaleIdx).numSuccess, ...
        statsByScale(scaleIdx).numTotal, ...
        statsByScale(scaleIdx).successRate, ...
        statsByScale(scaleIdx).numMetricValues, ...
        metricSpec.logName, median(values), numFloorClipped);
end

if foldYAxis
    format_folded_axes(topAx, bottomAx, scaleValues, lowerYLim, upperYLim, ...
        metricSpec.useLogScale, metricSpec.yLabel, 'Sampling region scale', ...
        'Sensitivity to initial data');
else
    set(ax, ...
        'XTick', 1:numel(scaleValues), ...
        'XTickLabel', compose('%.3g', scaleValues), ...
        'FontName', 'Arial', ...
        'FontSize', 11);
    if metricSpec.useLogScale
        set(ax, 'YScale', 'log');
    end
    xlim(ax, [0.5, numel(scaleValues) + 0.7]);
    ylim(ax, [displayFloor, 1.25 * max(max(allDistances), displayFloor)]);
    xlabel(ax, 'Sampling region scale', 'FontName', 'Arial');
    ylabel(ax, metricSpec.yLabel, 'FontName', 'Arial');
    title(ax, 'Sensitivity to initial data', ...
        'FontName', 'Arial', 'FontWeight', 'normal');
end

fprintf('Display floor=%.6g\n', displayFloor);

exportgraphics(fig, fullfile(resultsRootDir, metricSpec.outputFileName), ...
    'Resolution', 300);

function metricSpec = build_metric_spec(errorMetric, circuitDir)
switch lower(errorMetric)
    case 'parameter'
        baseData = load(fullfile(circuitDir, 'learnedData_ODE.mat'), 'Parameters');
        basePhysicalParams = to_physical_params(baseData.Parameters);
        distanceScale = abs(basePhysicalParams);
        distanceScale(distanceScale == 0) = 1;

        metricSpec = struct( ...
            'name', 'parameter', ...
            'statsFileName', 'successful_params_stats.mat', ...
            'outputFileName', 'Sensitivity_to_initial_data_violin.png', ...
            'yLabel', 'Relative parameter diff', ...
            'logName', 'RelativeDistance', ...
            'useLogScale', true, ...
            'useLogDensity', true, ...
            'basePhysicalParams', basePhysicalParams, ...
            'distanceScale', distanceScale);

    case 'orbit_l2'
        metricSpec = struct( ...
            'name', 'orbit_l2', ...
            'statsFileName', 'successful_params_stats.mat', ...
            'outputFileName', 'Sensitivity_to_initial_data_orbit_l2_violin.png', ...
            'yLabel', 'Relative L2 error', ...
            'logName', 'OrbitRelativeL2', ...
            'useLogScale', true, ...
            'useLogDensity', true, ...
            'basePhysicalParams', [], ...
            'distanceScale', []);
end
end

function values = collect_metric_values(statsData, metricSpec)
switch metricSpec.name
    case 'parameter'
        relativeDifferences = (statsData.successPhysicalParams - metricSpec.basePhysicalParams) ./ ...
            metricSpec.distanceScale;
        values = vecnorm(relativeDifferences, 2, 2);

    case 'orbit_l2'
        values = statsData.orbitRelativeL2(statsData.orbitSuccess);
end

values = values(:);
values = values(isfinite(values));
end

function displayFloor = choose_display_floor(values, useLogScale)
if useLogScale
    positiveValues = values(values > 0);
    displayFloor = 10^(floor(log10(prctile(positiveValues, 5))) - 1);
else
    displayFloor = min(values);
end
end

function [lowerYLim, upperYLim, lowerCut, upperCut] = folded_y_limits(values, padding, useLogScale)
[lowerCut, upperCut] = largest_gap_cut(values, useLogScale);
mainValues = values(values <= lowerCut);
outlierValues = values(values >= upperCut);

lowerYLim = padded_limits(mainValues, padding, useLogScale);
upperYLim = padded_limits(outlierValues, padding, useLogScale);
end

function limits = padded_limits(values, padding, useLogScale)
if useLogScale
    logValues = log10(values(values > 0));
    logMin = min(logValues);
    logMax = max(logValues);
    logSpan = max(logMax - logMin, 0.05);
    limits = 10 .^ [logMin - padding * logSpan, logMax + padding * logSpan];
else
    valueMin = min(values);
    valueMax = max(values);
    valueSpan = max(valueMax - valueMin, eps);
    limits = [valueMin - padding * valueSpan, valueMax + padding * valueSpan];
end
end

function [lowerCut, upperCut] = largest_gap_cut(values, useLogScale)
if useLogScale
    sortedValues = sort(log10(values(values > 0)));
    [~, gapIdx] = max(diff(sortedValues));
    lowerCut = 10 ^ sortedValues(gapIdx);
    upperCut = 10 ^ sortedValues(gapIdx + 1);
else
    sortedValues = sort(values);
    [~, gapIdx] = max(diff(sortedValues));
    lowerCut = sortedValues(gapIdx);
    upperCut = sortedValues(gapIdx + 1);
end
end

function foldYAxis = should_fold_y_axis(errorMetric)
foldYAxis = strcmpi(errorMetric, 'orbit_l2');
end

function physicalParams = to_physical_params(parameters)
parameters = reshape(double(parameters), 1, []);
physicalParams = [parameters(1), 1 / parameters(4), 1 / parameters(5), 1 / parameters(6)];
end

function draw_raincloud_group(ax, xCenter, values, color, useLogDensity)
values = values(:);
values = values(isfinite(values));

draw_scatter_group(ax, xCenter, values, color);
draw_half_violin(ax, xCenter + 0.04, values, color, useLogDensity);
draw_box_summary(ax, xCenter + 0.03, values, color);
end

function draw_scatter_group(ax, xCenter, values, color)
values = values(:);
values = values(isfinite(values));

scatterWidth = 0.16;
scatterX = xCenter - 0.26 + scatterWidth * (rand(size(values)) - 0.5);
scatter(ax, scatterX, values, 5, ...
    'MarkerFaceColor', color, ...
    'MarkerEdgeColor', 'w', ...
    'MarkerFaceAlpha', 0.65, ...
    'MarkerEdgeAlpha', 0.75, ...
    'LineWidth', 0.4);
end

function draw_half_violin(ax, xBase, values, color, useLogDensity)
if numel(unique(values)) < 2
    plot(ax, [xBase, xBase + 0.22], [values(1), values(1)], '-', ...
        'Color', lighten_color(color, 0.35), 'LineWidth', 6);
    return
end

if useLogDensity
    [density, logSupport] = ksdensity(log10(values(values > 0)), 'NumPoints', 150);
    support = 10 .^ logSupport;
else
    [density, support] = ksdensity(values, 'NumPoints', 150);
end

violinWidth = 0.28 * density / max(density);
xPatch = [xBase + violinWidth, repmat(xBase, 1, numel(support))];
yPatch = [support, fliplr(support)];
patch(ax, xPatch, yPatch, color, ...
    'FaceAlpha', 0.22, ...
    'EdgeColor', lighten_color(color, 0.45), ...
    'EdgeAlpha', 0.45, ...
    'LineWidth', 0.8);
end

function draw_box_summary(ax, xCenter, values, color)
quartiles = prctile(values, [25, 50, 75]);
q1 = quartiles(1);
medianValue = quartiles(2);
q3 = quartiles(3);
whiskerLow = min(values);
whiskerHigh = max(values);

boxWidth = 0.18;
boxLeft = xCenter - boxWidth / 2;
boxRight = xCenter + boxWidth / 2;
boxLineWidth = 0.5;
medianLineWidth = 0.8;

plot(ax, [xCenter, xCenter], [whiskerLow, q1], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', boxLineWidth);
plot(ax, [xCenter, xCenter], [q3, whiskerHigh], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', boxLineWidth);
plot(ax, [boxLeft, boxRight], [whiskerLow, whiskerLow], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', boxLineWidth);
plot(ax, [boxLeft, boxRight], [whiskerHigh, whiskerHigh], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', boxLineWidth);

rectangle(ax, ...
    'Position', [boxLeft, q1, boxWidth, max(q3 - q1, eps)], ...
    'FaceColor', [color, 0.58], ...
    'EdgeColor', 'none');
plot(ax, [boxLeft, boxRight], [medianValue, medianValue], 'k-', ...
    'LineWidth', medianLineWidth);
end

function lighterColor = lighten_color(color, amount)
lighterColor = color + amount * (1 - color);
lighterColor = min(max(lighterColor, 0), 1);
end

function [topAx, bottomAx] = create_folded_axes(fig)
bottomAx = axes(fig, 'Position', [0.12, 0.13, 0.82, 0.55]);
hold(bottomAx, 'on')
grid(bottomAx, 'on')
box(bottomAx, 'on')

topAx = axes(fig, 'Position', [0.12, 0.74, 0.82, 0.18]);
hold(topAx, 'on')
grid(topAx, 'on')
box(topAx, 'on')
end

function format_folded_axes(topAx, bottomAx, xValues, lowerYLim, upperYLim, ...
    useLogScale, yLabelText, xLabelText, titleText)
axesList = [topAx, bottomAx];
for axIdx = 1:numel(axesList)
    ax = axesList(axIdx);
    set(ax, ...
        'XTick', 1:numel(xValues), ...
        'FontName', 'Arial', ...
        'FontSize', 11);
    if useLogScale
        set(ax, 'YScale', 'log');
    end
    xlim(ax, [0.5, numel(xValues) + 0.7]);
end

ylim(bottomAx, lowerYLim);
ylim(topAx, upperYLim);
set(bottomAx, 'XTickLabel', compose('%.3g', xValues));
set(topAx, 'XTickLabel', []);
xlabel(bottomAx, xLabelText, 'FontName', 'Arial');
ylabel(bottomAx, yLabelText, 'FontName', 'Arial');
title(topAx, titleText, 'FontName', 'Arial', 'FontWeight', 'normal');

text(bottomAx, -0.055, 1.02, '//', 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
text(topAx, -0.055, -0.10, '//', 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end
