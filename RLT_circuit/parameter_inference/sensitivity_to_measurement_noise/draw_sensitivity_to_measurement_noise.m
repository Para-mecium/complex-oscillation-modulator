clear
clc

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);

baseFile = fullfile(circuitDir, 'learnedData_ODE.mat');
resultsDir = fullfile(scriptDir, 'modulation_results');

if ~isfile(baseFile)
    error('draw_sensitivity_to_measurement_noise:MissingBaseline', ...
        'Baseline learned data not found: %s. Run params_inf.m first.', baseFile);
end
if ~isfolder(resultsDir)
    error('draw_sensitivity_to_measurement_noise:MissingResults', ...
        'Modulation results directory not found: %s.', resultsDir);
end

baseData = load(baseFile, 'Parameters');
basePhysicalParams = to_physical_params(baseData.Parameters);
distanceScale = abs(basePhysicalParams);
distanceScale(distanceScale == 0) = 1;

summaryFiles = dir(fullfile(resultsDir, 'noise_level_*', 'modulation_summary.mat'));
if isempty(summaryFiles)
    error('draw_sensitivity_to_measurement_noise:NoSummaryFiles', ...
        'No modulation_summary.mat files found under %s.', resultsDir);
end

noiseLevels = zeros(numel(summaryFiles), 1);
distancesByLevel = cell(numel(summaryFiles), 1);
successRates = zeros(numel(summaryFiles), 1);

for fileIdx = 1:numel(summaryFiles)
    summaryFile = fullfile(summaryFiles(fileIdx).folder, summaryFiles(fileIdx).name);
    summaryData = load(summaryFile, 'noiseLevel', 'sampleSuccess', 'identifiedParameters', 'successRate');

    noiseLevels(fileIdx) = summaryData.noiseLevel;
    successRates(fileIdx) = summaryData.successRate;

    sampleSuccess = reshape(summaryData.sampleSuccess, 1, []);
    identifiedParameters = reshape(summaryData.identifiedParameters, 1, []);

    distances = [];
    for sampleIdx = 1:numel(identifiedParameters)
        if sampleIdx > numel(sampleSuccess) || ~sampleSuccess(sampleIdx)
            continue
        end

        parameters = identifiedParameters{sampleIdx};
        if isempty(parameters) || any(~isfinite(parameters(:)))
            continue
        end

        physicalParams = to_physical_params(parameters);
        relativeDifference = (physicalParams - basePhysicalParams) ./ distanceScale;
        distances(end + 1, 1) = norm(relativeDifference, 2); %#ok<SAGROW>
    end

    distancesByLevel{fileIdx} = distances;
end

[noiseLevels, sortIdx] = sort(noiseLevels);
distancesByLevel = distancesByLevel(sortIdx);
successRates = successRates(sortIdx);

fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on')
grid(ax, 'on')
box(ax, 'on')

colors = lines(numel(noiseLevels));
rng(1, 'twister');

for levelIdx = 1:numel(noiseLevels)
    values = distancesByLevel{levelIdx};
    if isempty(values)
        warning('No successful samples for noiseLevel=%.6g.', noiseLevels(levelIdx));
        continue
    end

    draw_raincloud_group(ax, levelIdx, values, colors(levelIdx, :));
    fprintf('noiseLevel=%.6g, successRate=%.3f, n=%d, medianRelativeDistance=%.6g\n', ...
        noiseLevels(levelIdx), successRates(levelIdx), numel(values), median(values));
end

set(ax, ...
    'XTick', 1:numel(noiseLevels), ...
    'XTickLabel', compose('%.3g', noiseLevels), ...
    'FontName', 'Arial', ...
    'FontSize', 11);
xlim(ax, [0.5, numel(noiseLevels) + 0.7]);
xlabel(ax, 'Noise intensity', 'FontName', 'Arial');
ylabel(ax, 'Relative Euclidean distance from noiseless inferred physical parameters', ...
    'FontName', 'Arial');
title(ax, 'Sensitivity to measurement noise', 'FontName', 'Arial', 'FontWeight', 'normal');

function physicalParams = to_physical_params(parameters)
parameters = reshape(double(parameters), 1, []);
if numel(parameters) < 6
    error('draw_sensitivity_to_measurement_noise:InvalidParameterVector', ...
        'Expected at least 6 parameters, received %d.', numel(parameters));
end

physicalParams = [parameters(1), 1 / parameters(4), 1 / parameters(5), 1 / parameters(6)];
end

function draw_raincloud_group(ax, xCenter, values, color)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    return
end

scatterWidth = 0.16;
scatterX = xCenter - 0.26 + scatterWidth * (rand(size(values)) - 0.5);
scatter(ax, scatterX, values, 28, ...
    'MarkerFaceColor', color, ...
    'MarkerEdgeColor', 'w', ...
    'MarkerFaceAlpha', 0.65, ...
    'MarkerEdgeAlpha', 0.75, ...
    'LineWidth', 0.4);

    draw_half_violin(ax, xCenter + 0.04, values, color);
    draw_box_summary(ax, xCenter + 0.03, values, color);
end

function draw_half_violin(ax, xBase, values, color)
if numel(unique(values)) < 2
    plot(ax, [xBase, xBase + 0.22], [values(1), values(1)], '-', ...
        'Color', lighten_color(color, 0.35), 'LineWidth', 6);
    return
end

[density, support] = ksdensity(values, 'NumPoints', 150);
if all(density == 0)
    return
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

plot(ax, [xCenter, xCenter], [whiskerLow, q1], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 1.8);
plot(ax, [xCenter, xCenter], [q3, whiskerHigh], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 1.8);
plot(ax, [boxLeft, boxRight], [whiskerLow, whiskerLow], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 1.8);
plot(ax, [boxLeft, boxRight], [whiskerHigh, whiskerHigh], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 1.8);

rectangle(ax, ...
    'Position', [boxLeft, q1, boxWidth, max(q3 - q1, eps)], ...
    'FaceColor', [color, 0.58], ...
    'EdgeColor', 'none');
plot(ax, [boxLeft, boxRight], [medianValue, medianValue], 'k-', 'LineWidth', 1.8);
end

function lighterColor = lighten_color(color, amount)
lighterColor = color + amount * (1 - color);
lighterColor = min(max(lighterColor, 0), 1);
end
