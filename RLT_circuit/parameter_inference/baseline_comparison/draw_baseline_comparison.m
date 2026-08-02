clear
clc
%%
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
resultsDir = fullfile(scriptDir, 'results');
lossName = 'relative_l2_orbit';
proposedSensitivityScale = 1;
fontSize = 9;
numTargetScatterPoints = 50;
proposedLineWidth = 0.5;
baselineLineWidth = 0.5;
baselineLineStyles = {'-', '-', '--', '-.', ':'};
addpath(fullfile(repoDir, 'PO_extract'), '-begin');
summaryFile = fullfile(resultsDir, ...
    sprintf('baseline_comparison_summary_%s.mat', result_name_token(lossName)));
legacySummaryFile = fullfile(resultsDir, 'baseline_comparison_summary.mat');
if ~isfile(summaryFile) && isfile(legacySummaryFile)
    summaryFile = legacySummaryFile;
end

if ~isfile(summaryFile)
    error('draw_baseline_comparison:MissingSummary', ...
        'Run run_baseline_comparison.m first. Missing file: %s', summaryFile);
end

summaryData = load(summaryFile);
summaryRows = summaryData.summaryRows;
config = build_config(summaryData.configSummary);
lossName = config.lossName;
targetOrbit = config.targetOrbit;

records = [ ...
    load_proposed_records(resultsDir, lossName), ...
    load_result_records(summaryRows, resultsDir, lossName)];
if isempty(records)
    error('draw_baseline_comparison:NoResults', ...
        'No successful baseline result files found in %s.', resultsDir);
end

[groupKeys, groupLabels, groupIndex] = unique_groups(records);
colors = lines(numel(groupKeys));
sobolGroup = startsWith(groupKeys, 'sobol|');
if any(sobolGroup)
    colors(sobolGroup, :) = repmat([0.2000, 0.2000, 0.2000], nnz(sobolGroup), 1);
end
proposedGroup = is_proposed_group(groupKeys, groupLabels);
axisPadding = 0.08;
finalLossesByGroup = group_record_values(records, groupIndex, numel(groupKeys), 'finalLoss');
runtimesByGroup = group_record_values(records, groupIndex, numel(groupKeys), 'runtime');
figure2FinalLossesByGroup = finalLossesByGroup;
runtimeFinalLossesByGroup = finalLossesByGroup;
runtimeRuntimesByGroup = runtimesByGroup;
proposedSensitivity = load_proposed_sensitivity_runtime_loss( ...
    parameterInferenceDir, config, proposedSensitivityScale);
proposedGroupIdx = find(proposedGroup, 1, 'first');
if ~isempty(proposedGroupIdx) && ~isempty(proposedSensitivity.successRuntime)
    figure2FinalLossesByGroup{proposedGroupIdx} = proposedSensitivity.successLoss;
    runtimeFinalLossesByGroup{proposedGroupIdx} = proposedSensitivity.successLoss;
    runtimeRuntimesByGroup{proposedGroupIdx} = proposedSensitivity.successRuntime;
end

%% Figure 1: Best-so-far loss vs evaluations
figure
hold on
grid on
box on

for groupIdx = 1:numel(groupKeys)
    if proposedGroup(groupIdx)
        continue
    end

    groupRecords = records(groupIndex == groupIdx);
    curves = collect_curves(groupRecords);
    evalAxis = 1:size(curves, 2);
    medianCurve = row_nan_percentile(curves, 50);
    q1Curve = row_nan_percentile(curves, 25);
    q3Curve = row_nan_percentile(curves, 75);

    semilogy(evalAxis, medianCurve, '-', 'LineWidth', 1, ...
        'Color', colors(groupIdx, :), 'DisplayName', groupLabels{groupIdx});
    semilogy(evalAxis, q1Curve, '--', 'LineWidth', 0.5, ...
        'Color', colors(groupIdx, :), 'HandleVisibility', 'off');
    semilogy(evalAxis, q3Curve, '--', 'LineWidth', 0.5, ...
        'Color', colors(groupIdx, :), 'HandleVisibility', 'off');
end

xlabel('Forward evaluations')
ylabel('Best-so-far loss')
set(gca, 'FontName', 'Arial', 'FontSize', fontSize, 'YScale', 'log')
ylim(plotted_y_limits(gca, axisPadding, true))

%% Figure 2: Final best loss across seeds
fig2 = figure('Color', 'w');
ax2 = axes(fig2);
hold(ax2, 'on')
grid(ax2, 'on')
box(ax2, 'on')
rng(1, 'twister')
for groupIdx = 1:numel(groupKeys)
    draw_distribution_violin(ax2, groupIdx, ...
        figure2FinalLossesByGroup{groupIdx}, colors(groupIdx, :));
end
set(ax2, ...
    'XTick', 1:numel(groupLabels), ...
    'FontName', 'Arial', ...
    'FontSize', fontSize, ...
    'YScale', 'log')
xlim(ax2, [0.5, numel(groupLabels) + 0.5])
ylim(ax2, plotted_y_limits(ax2, axisPadding, true))
ylabel(ax2, 'Final best loss', 'FontName', 'Arial')

%% Figure 3: Runtime vs final best loss
fig3 = figure('Color', 'w');
[mainAx3, failureAx3] = create_runtime_axes(fig3);
for groupIdx = 1:numel(groupKeys)
    finalLosses = runtimeFinalLossesByGroup{groupIdx};
    runtimes = runtimeRuntimesByGroup{groupIdx};
    mask = isfinite(finalLosses) & finalLosses > 0 & ...
        isfinite(runtimes) & runtimes > 0;
    scatter(mainAx3, runtimes(mask), finalLosses(mask), 5, colors(groupIdx, :), ...
        'filled', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.65, ...
        'HandleVisibility', 'off');
end
if ~isempty(proposedGroupIdx) && ~isempty(proposedSensitivity.failedRuntime)
    scatter(failureAx3, proposedSensitivity.failedRuntime, ...
        0.5 * ones(size(proposedSensitivity.failedRuntime)), 5, ...
        colors(proposedGroupIdx, :), 'x', 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end
runtimeXValues = [vertcat(runtimeRuntimesByGroup{:}); proposedSensitivity.failedRuntime(:)];
runtimeYValues = vertcat(runtimeFinalLossesByGroup{:});
runtimeYValues = runtimeYValues(isfinite(runtimeYValues) & runtimeYValues > 0);
format_runtime_axes(mainAx3, failureAx3, runtimeXValues, runtimeYValues, axisPadding, fontSize, ...
    'Runtime (s)', 'Final best loss');

%% Figure 4: Best parameter time series comparison
figure
stateCount = size(targetOrbit.y, 2);
representatives = select_group_representatives(records, groupIndex, numel(groupKeys));
plotPhase = linspace(0, 1, config.lossOptions.compareNumPoints + 1).';
plotPhase(end) = [];
targetY = resample_orbit_y(targetOrbit, plotPhase);
targetScatterPhase = linspace(0, 1, numTargetScatterPoints + 1).';
targetScatterPhase(end) = [];
targetScatterY = resample_orbit_y(targetOrbit, targetScatterPhase);

for stateIdx = 1:stateCount
    subplot(stateCount, 1, stateIdx)
    hold on
    grid on
    box on

    baselineStyleIdx = 0;
    for groupIdx = 1:numel(groupKeys)
        record = representatives(groupIdx);
        if isempty(record.bestOrbit)
            continue
        end

        candidateY = resample_orbit_y(record.bestOrbit, plotPhase);
        candidateY = align_candidate_phase(candidateY, targetY);
        if proposedGroup(groupIdx)
            lineColor = colors(groupIdx, :);
            lineWidth = proposedLineWidth;
            lineStyle = '-';
        else
            baselineStyleIdx = baselineStyleIdx + 1;
            lineColor = colors(groupIdx, :);
            lineWidth = baselineLineWidth;
            lineStyle = baselineLineStyles{baselineStyleIdx};
        end
        plot(plotPhase, candidateY(:, stateIdx), 'LineWidth', lineWidth, ...
            'Color', lineColor, ...
            'LineStyle', lineStyle, ...
            'DisplayName', groupLabels{groupIdx});
    end

    scatter(targetScatterPhase, targetScatterY(:, stateIdx), 5, 'k', 'x', ...
        'LineWidth', 0.6, 'DisplayName', 'target','MarkerEdgeColor','r');

    ylabel(sprintf('V_{%d}', stateIdx))
    if stateIdx == 1
        ylim([0 2])
    end
    if stateIdx == 0
        lgd = legend(gca, 'Location', 'northeast', ...
            'NumColumns', 3, 'Box', 'off');
        set(lgd, 'FontName', 'Arial', 'FontSize', 6, ...
            'ItemTokenSize', [8, 4]);
    end
    if stateIdx == stateCount
        xlabel('Normalized phase')
    end
    set(gca, 'FontName', 'Arial', 'FontSize', fontSize)
end

%% Figure 5: Best parameter residual time series comparison
figure
for stateIdx = 1:stateCount
    subplot(stateCount, 1, stateIdx)
    hold on
    grid on
    box on
    yline(0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');

    baselineStyleIdx = 0;
    for groupIdx = 1:numel(groupKeys)
        record = representatives(groupIdx);

        candidateY = resample_orbit_y(record.bestOrbit, plotPhase);
        candidateY = align_candidate_phase(candidateY, targetY);
        residualY = candidateY - targetY;
        if proposedGroup(groupIdx)
            lineColor = colors(groupIdx, :);
            lineWidth = proposedLineWidth;
            lineStyle = '-';
        else
            baselineStyleIdx = baselineStyleIdx + 1;
            lineColor = colors(groupIdx, :);
            lineWidth = baselineLineWidth;
            lineStyle = baselineLineStyles{baselineStyleIdx};
        end
        plot(plotPhase, residualY(:, stateIdx), 'LineWidth', lineWidth, ...
            'Color', lineColor, ...
            'LineStyle', lineStyle, ...
            'DisplayName', groupLabels{groupIdx});
    end

    ylabel(sprintf('\\Delta V_{%d}', stateIdx))
    if stateIdx == 1
        lgd = legend(gca, 'Location', 'northeast', ...
            'NumColumns', 3, 'Box', 'off');
        set(lgd, 'FontName', 'Arial', 'FontSize', 6, ...
            'ItemTokenSize', [8, 4]);
    end
    if stateIdx == stateCount
        xlabel('Normalized phase')
    end
    set(gca, 'FontName', 'Arial', 'FontSize', fontSize)
end

%% Local functions
function records = load_result_records(summaryRows, resultsDir, lossName)
records = struct( ...
    'methodName', {}, ...
    'refinementMethod', {}, ...
    'groupKey', {}, ...
    'groupLabel', {}, ...
    'seed', {}, ...
    'finalLoss', {}, ...
    'runtime', {}, ...
    'bestSoFarLosses', {}, ...
    'bestActiveParams', {}, ...
    'bestOrbit', {});

for i = 1:numel(summaryRows)
    row = summaryRows(i);
    if ~strcmp(row.status, 'success') || isempty(row.resultFile)
        continue
    end

    resultFile = resolve_result_file(row, resultsDir, lossName);
    if ~isfile(resultFile)
        continue
    end

    fileData = load(resultFile);
    result = fileData.result;
    if isfield(result, 'refinement') && isfield(result.refinement, 'method')
        refinementMethod = result.refinement.method;
    else
        refinementMethod = 'none';
    end
    groupLabel = make_group_label(result.methodName, refinementMethod);

    record = struct();
    record.methodName = result.methodName;
    record.refinementMethod = refinementMethod;
    record.groupKey = sprintf('%s|%s', result.methodName, refinementMethod);
    record.groupLabel = groupLabel;
    record.seed = result.seed;
    record.finalLoss = final_loss(result);
    record.runtime = result_runtime(result);
    record.bestSoFarLosses = best_so_far_curve(result);
    [record.bestActiveParams, record.bestOrbit] = final_best_candidate(result);

    records(end + 1) = record; %#ok<AGROW>
end
end

function resultFile = resolve_result_file(row, resultsDir, lossName)
resultFile = row.resultFile;
if isfile(resultFile)
    return
end

[~, fileName, fileExt] = fileparts(resultFile);
candidateFile = fullfile( ...
    result_method_loss_dir(resultsDir, row.methodName, lossName), ...
    [fileName fileExt]);
if isfile(candidateFile)
    resultFile = candidateFile;
end
end

function records = load_proposed_records(resultsDir, lossName)
records = struct( ...
    'methodName', {}, ...
    'refinementMethod', {}, ...
    'groupKey', {}, ...
    'groupLabel', {}, ...
    'seed', {}, ...
    'finalLoss', {}, ...
    'runtime', {}, ...
    'bestSoFarLosses', {}, ...
    'bestActiveParams', {}, ...
    'bestOrbit', {});

proposedFile = fullfile( ...
    result_method_loss_dir(resultsDir, 'proposed_method', lossName), ...
    'proposed_method_result.mat');
if ~isfile(proposedFile)
    return
end

fileData = load(proposedFile);
result = fileData.result;

record = struct();
record.methodName = result.methodName;
record.refinementMethod = 'none';
record.groupKey = sprintf('%s|none', result.methodName);
record.groupLabel = 'proposed method';
record.seed = NaN;
record.finalLoss = final_loss(result);
record.runtime = result_runtime(result);
record.bestSoFarLosses = best_so_far_curve(result);
[record.bestActiveParams, record.bestOrbit] = final_best_candidate(result);

records = record;
end

function label = make_group_label(methodName, refinementMethod)
if strcmp(refinementMethod, 'none')
    label = methodName;
else
    label = sprintf('%s+%s', methodName, refinementMethod);
end
end

function loss = final_loss(result)
if isfield(result, 'refinement') && result.refinement.enabled
    loss = result.refinement.bestLoss;
else
    loss = result.bestLoss;
end
end

function runtime = result_runtime(result)
if isfield(result, 'runtime') && isfinite(result.runtime)
    runtime = result.runtime;
elseif isfield(result, 'searchRuntime') && isfield(result, 'refinementRuntime')
    runtime = result.searchRuntime + result.refinementRuntime;
else
    runtime = NaN;
end
end

function curve = best_so_far_curve(result)
if isfield(result, 'totalBestSoFarLosses') && ~isempty(result.totalBestSoFarLosses)
    curve = result.totalBestSoFarLosses(:).';
    return
end

if isfield(result, 'bestSoFarLosses') && ~isempty(result.bestSoFarLosses)
    searchCurve = result.bestSoFarLosses(:);
else
    searchCurve = cummin(result.losses(:));
end

if isfield(result, 'refinement') && result.refinement.enabled && ...
        isfield(result.refinement, 'globalBestSoFarLosses')
    curve = [searchCurve; result.refinement.globalBestSoFarLosses(:)].';
else
    curve = searchCurve.';
end
end

function [bestActiveParams, bestOrbit] = final_best_candidate(result)
bestActiveParams = result.bestActiveParams;
bestOrbit = [];

if isfield(result, 'refinement') && result.refinement.enabled && ...
        result.refinement.bestLoss <= result.bestLoss
    bestActiveParams = result.refinement.bestActiveParams;
    if ~isempty(result.refinement.bestForwardResult) && ...
            isfield(result.refinement.bestForwardResult, 'orbit')
        bestOrbit = result.refinement.bestForwardResult.orbit;
    end
elseif ~isempty(result.bestForwardResult) && isfield(result.bestForwardResult, 'orbit')
    bestOrbit = result.bestForwardResult.orbit;
end
end

function [groupKeys, groupLabels, groupIndex] = unique_groups(records)
allKeys = {records.groupKey};
[groupKeys, ~, groupIndex] = unique(allKeys, 'stable');
groupLabels = cell(size(groupKeys));
for i = 1:numel(groupKeys)
    firstIdx = find(groupIndex == i, 1, 'first');
    groupLabels{i} = records(firstIdx).groupLabel;
end
end

function curves = collect_curves(groupRecords)
maxLength = max(arrayfun(@(r) numel(r.bestSoFarLosses), groupRecords));
curves = NaN(numel(groupRecords), maxLength);
for i = 1:numel(groupRecords)
    curve = groupRecords(i).bestSoFarLosses(:).';
    curves(i, 1:numel(curve)) = curve;
end
end

function values = row_nan_percentile(curves, percentile)
values = NaN(1, size(curves, 2));
for j = 1:size(curves, 2)
    column = curves(:, j);
    column = column(isfinite(column));
    if isempty(column)
        continue
    end
    values(j) = percentile_value(column, percentile);
end
end

function value = percentile_value(data, percentile)
data = sort(data(:));
if isscalar(data)
    value = data;
    return
end
position = 1 + (numel(data) - 1) * percentile / 100;
lowerIdx = floor(position);
upperIdx = ceil(position);
if lowerIdx == upperIdx
    value = data(lowerIdx);
else
    weight = position - lowerIdx;
    value = (1 - weight) * data(lowerIdx) + weight * data(upperIdx);
end
end

function proposedGroup = is_proposed_group(groupKeys, groupLabels)
proposedGroup = false(size(groupKeys));
for groupIdx = 1:numel(groupKeys)
    proposedGroup(groupIdx) = strcmp(groupLabels{groupIdx}, 'proposed method') || ...
        startsWith(groupKeys{groupIdx}, 'proposed_method|');
end
end

function valuesByGroup = group_record_values(records, groupIndex, numGroups, fieldName)
valuesByGroup = cell(numGroups, 1);
for groupIdx = 1:numGroups
    groupRecords = records(groupIndex == groupIdx);
    valuesByGroup{groupIdx} = reshape([groupRecords.(fieldName)], [], 1);
end
end

function proposedData = load_proposed_sensitivity_runtime_loss( ...
    parameterInferenceDir, config, scale)
sensitivityDir = fullfile(parameterInferenceDir, 'sensitivity_to_init_data');
scaleDir = fullfile(sensitivityDir, 'results', scale_dir_name(scale));
statsFile = fullfile(scaleDir, 'successful_params_stats.mat');
summaryFile = fullfile(scaleDir, 'params_inf_sensitivity_summary.mat');

statsData = load(statsFile);
summaryData = load(summaryFile);

numSuccess = numel(statsData.successFinalOrbit);
successLoss = NaN(numSuccess, 1);
for sampleIdx = 1:numSuccess
    if statsData.orbitSuccess(sampleIdx)
        candidateOrbit = statsData.successFinalOrbit{sampleIdx};
        candidateFeatures = evaluate_orbit_features(candidateOrbit, [], [], struct());
        successLoss(sampleIdx) = loss_function(candidateOrbit, ...
            config.targetOrbit, config.lossOptions, ...
            candidateFeatures, config.targetFeatures);
    end
end

successRuntime = statsData.successRuntimeFMAM(:);
successMask = isfinite(successRuntime) & successRuntime > 0 & ...
    isfinite(successLoss) & successLoss > 0;

runResults = summaryData.runResults;
failedMask = ~[runResults.success];
failedRuntime = [runResults.runtimeFMAM].';
proposedData = struct();
proposedData.successRuntime = successRuntime(successMask);
proposedData.successLoss = successLoss(successMask);
proposedData.failedRuntime = failedRuntime(failedMask(:) & isfinite(failedRuntime) & failedRuntime > 0);
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end

function [mainAx, failureAx] = create_runtime_axes(fig)
failureAx = axes(fig, 'Position', [0.12, 0.15, 0.82, 0.08]);
hold(failureAx, 'on')
grid(failureAx, 'on')
box(failureAx, 'on')

mainAx = axes(fig, 'Position', [0.12, 0.25, 0.82, 0.68]);
hold(mainAx, 'on')
grid(mainAx, 'on')
box(mainAx, 'on')
end

function format_runtime_axes(mainAx, failureAx, xValues, yValues, axisPadding, fontSize, ...
    xLabelText, yLabelText)
xLimits = padded_log_limits(xValues, axisPadding);
set(mainAx, ...
    'FontName', 'Arial', ...
    'FontSize', fontSize, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'XTickLabel', [])
set(failureAx, ...
    'FontName', 'Arial', ...
    'FontSize', fontSize, ...
    'XScale', 'log', ...
    'YTick', 0.5, ...
    'YTickLabel', {'failed'})
xlim(mainAx, xLimits)
xlim(failureAx, xLimits)
ylim(mainAx, padded_limits(yValues, axisPadding, true))
ylim(failureAx, [0 1])
ylabel(mainAx, yLabelText, 'FontName', 'Arial')
xlabel(failureAx, xLabelText, 'FontName', 'Arial')
end

function limits = plotted_y_limits(ax, padding, useLogScale)
children = findobj(ax, '-property', 'YData');
values = [];
for childIdx = 1:numel(children)
    yData = get(children(childIdx), 'YData');
    values = [values; yData(:)]; %#ok<AGROW>
end

values = values(isfinite(values));
if useLogScale
    values = values(values > 0);
end
limits = padded_limits(values, padding, useLogScale);
end

function limits = padded_limits(values, padding, useLogScale)
values = values(:);
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

function limits = padded_log_limits(values, padding)
values = values(:);
values = values(isfinite(values) & values > 0);
logValues = log10(values);
logMin = min(logValues);
logMax = max(logValues);
logSpan = max(logMax - logMin, 0.05);
limits = 10 .^ [logMin - padding * logSpan, logMax + padding * logSpan];
end

function draw_distribution_violin(ax, x, values, color)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    return
end

draw_distribution_scatter(ax, x, values, color);
if numel(unique(values)) >= 2
    draw_half_violin(ax, x + 0.04, values, color);
    draw_box_summary(ax, x + 0.03, values, color);
end
end

function draw_distribution_scatter(ax, xCenter, values, color)
scatterWidth = 0.16;
if isscalar(values)
    scatterX = xCenter;
else
    scatterX = xCenter - 0.26 + scatterWidth * (rand(size(values)) - 0.5);
end
scatter(ax, scatterX, values, 5, color, 'filled', ...
     'MarkerEdgeColor', 'w', ...
    'MarkerFaceAlpha', 0.65, ...
    'MarkerEdgeAlpha', 0.75, ...
    'LineWidth', 0.4);
end

function draw_half_violin(ax, xBase, values, color)
[density, logSupport] = ksdensity(log10(values(values > 0)), 'NumPoints', 150);
support = 10 .^ logSupport;
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
q1 = percentile_value(values, 25);
q2 = percentile_value(values, 50);
q3 = percentile_value(values, 75);
vMin = min(values);
vMax = max(values);
boxWidth = 0.18;
boxLeft = xCenter - boxWidth / 2;
boxRight = xCenter + boxWidth / 2;

plot(ax, [xCenter, xCenter], [vMin q1], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 0.5);
plot(ax, [xCenter, xCenter], [q3 vMax], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 0.5);
plot(ax, [boxLeft, boxRight], [vMin vMin], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 0.5);
plot(ax, [boxLeft, boxRight], [vMax vMax], '-', ...
    'Color', lighten_color(color, 0.2), 'LineWidth', 0.5);
rectangle(ax, ...
    'Position', [boxLeft, q1, boxWidth, max(q3 - q1, eps)], ...
    'FaceColor', [color, 0.58], ...
    'EdgeColor', 'none');
plot(ax, [boxLeft, boxRight], [q2 q2], 'k-', 'LineWidth', 0.8);
end

function lighterColor = lighten_color(color, amount)
lighterColor = color + amount * (1 - color);
lighterColor = min(max(lighterColor, 0), 1);
end

function representatives = select_group_representatives(records, groupIndex, numGroups)
representatives = records(ones(1, numGroups));
for groupIdx = 1:numGroups
    idx = find(groupIndex == groupIdx);
    losses = [records(idx).finalLoss];
    [~, bestLocalIdx] = min(losses);
    representatives(groupIdx) = records(idx(bestLocalIdx));
end
end

function yResampled = resample_orbit_y(orbit, phase)
orbitPhase = normalized_phase(orbit.t);
yResampled = interp1(orbitPhase, orbit.y, phase, 'linear');
end

function phase = normalized_phase(t)
t = t(:);
if numel(t) < 2 || t(end) == t(1)
    error('draw_baseline_comparison:InvalidOrbitTime', ...
        'Orbit time vector must contain at least two distinct values.');
end
phase = (t - t(1)) / (t(end) - t(1));
end

function alignedCandidateY = align_candidate_phase(candidateY, targetY)
nPoints = size(candidateY, 1);
targetNorm = norm(targetY(:));
if targetNorm == 0
    targetNorm = 1;
end

bestShift = 0;
bestLoss = Inf;
for shiftIdx = 0:(nPoints - 1)
    shiftedY = circshift(candidateY, shiftIdx, 1);
    shiftLoss = norm(shiftedY(:) - targetY(:)) / targetNorm;
    if shiftLoss < bestLoss
        bestLoss = shiftLoss;
        bestShift = shiftIdx;
    end
end

alignedCandidateY = circshift(candidateY, bestShift, 1);
end
