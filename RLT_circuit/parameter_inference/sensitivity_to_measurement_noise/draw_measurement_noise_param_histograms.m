clear
clc

scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'modulation_results');

if ~isfolder(resultsDir)
    error('draw_measurement_noise_param_histograms:MissingResults', ...
        'Modulation results directory not found: %s.', resultsDir);
end

summaryFiles = dir(fullfile(resultsDir, 'noise_level_*', 'modulation_summary.mat'));
if isempty(summaryFiles)
    error('draw_measurement_noise_param_histograms:NoSummaryFiles', ...
        'No modulation_summary.mat files found under %s.', resultsDir);
end

noiseLevels = zeros(numel(summaryFiles), 1);
successRates = zeros(numel(summaryFiles), 1);
physicalParamsByLevel = cell(numel(summaryFiles), 1);
noiseDirNames = cell(numel(summaryFiles), 1);

for fileIdx = 1:numel(summaryFiles)
    summaryFile = fullfile(summaryFiles(fileIdx).folder, summaryFiles(fileIdx).name);
    summaryData = load(summaryFile, ...
        'noiseLevel', 'sampleSuccess', 'identifiedParameters', 'successRate');

    noiseLevels(fileIdx) = summaryData.noiseLevel;
    successRates(fileIdx) = summaryData.successRate;
    [~, noiseDirNames{fileIdx}] = fileparts(summaryFiles(fileIdx).folder);
    physicalParamsByLevel{fileIdx} = collect_successful_physical_params( ...
        summaryData.identifiedParameters, summaryData.sampleSuccess);
end

[noiseLevels, sortIdx] = sort(noiseLevels);
successRates = successRates(sortIdx);
physicalParamsByLevel = physicalParamsByLevel(sortIdx);
noiseDirNames = noiseDirNames(sortIdx);

for levelIdx = 1:numel(noiseLevels)
    physicalParams = physicalParamsByLevel{levelIdx};
    if isempty(physicalParams)
        warning('No successful samples for noiseLevel=%.6g.', noiseLevels(levelIdx));
        continue
    end

    fprintf('noiseLevel=%.6g, successRate=%.3f, n=%d\n', ...
        noiseLevels(levelIdx), successRates(levelIdx), size(physicalParams, 1));
    fprintf('mean:   R_C=%.6g, C_1=%.6g, C_2=%.6g, C_3=%.6g\n', ...
        mean(physicalParams, 1));
    fprintf('median: R_C=%.6g, C_1=%.6g, C_2=%.6g, C_3=%.6g\n', ...
        median(physicalParams, 1));

    fig = draw_histogram_panel(physicalParams, noiseLevels(levelIdx));
    saveName = sprintf('successful_params_histograms_%s.png', noiseDirNames{levelIdx});
    exportgraphics(fig, fullfile(scriptDir, saveName), 'Resolution', 300);
    close(fig);
end

function physicalParams = collect_successful_physical_params(identifiedParameters, sampleSuccess)
sampleSuccess = reshape(sampleSuccess, 1, []);
identifiedParameters = reshape(identifiedParameters, 1, []);

physicalParams = zeros(0, 4);
for sampleIdx = 1:numel(identifiedParameters)
    if sampleIdx > numel(sampleSuccess) || ~sampleSuccess(sampleIdx)
        continue
    end

    parameters = identifiedParameters{sampleIdx};
    if isempty(parameters) || any(~isfinite(parameters(:)))
        continue
    end

    currentPhysicalParams = to_physical_params(parameters);
    if any(~isfinite(currentPhysicalParams))
        continue
    end

    physicalParams(end + 1, :) = currentPhysicalParams; %#ok<AGROW>
end
end

function fig = draw_histogram_panel(physicalParams, ~)
paramNames = {'R_C', 'C_1', 'C_2', 'C_3'};
xLabels = {'{\itR_C} (k\Omega)', ...
    '{\itC}_1 (\muF)', ...
    '{\itC}_2 (\muF)', ...
    '{\itC}_3 (\muF)'};

fig = figure('Color', 'w', 'Position', [100, 100, 900, 700]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for paramIdx = 1:4
    ax = nexttile;
    values = physicalParams(:, paramIdx);
    histogram(ax, values, 10, ...
        'FaceColor', [0.72 0.80 0.90], ...
        'EdgeColor', [0.25 0.25 0.25]);
    hold(ax, 'on')
    grid(ax, 'on')
    box(ax, 'on')

    xline(ax, mean(values), '-', ...
        'LineWidth', 1.2, ...
        'Color', [0.85 0.33 0.10]);
    xline(ax, median(values), '--', ...
        'LineWidth', 1.2, ...
        'Color', [0.00 0.45 0.74]);

    xlabel(ax, xLabels{paramIdx}, 'FontName', 'Arial')
    ylabel(ax, 'Count', 'FontName', 'Arial')
    title(ax, paramNames{paramIdx}, 'FontName', 'Arial')
    set(ax, 'FontName', 'Arial', 'FontSize', 10)
end
end

function physicalParams = to_physical_params(parameters)
parameters = reshape(double(parameters), 1, []);
if numel(parameters) < 6
    error('draw_measurement_noise_param_histograms:InvalidParameterVector', ...
        'Expected at least 6 parameters, received %d.', numel(parameters));
end

physicalParams = [parameters(1), 1 / parameters(4), 1 / parameters(5), 1 / parameters(6)];
end
