clear
clc

scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'results');
summaryFile = fullfile(resultsDir, 'baseline_comparison_summary.mat');

if ~isfile(summaryFile)
    error('draw_baseline_comparison:MissingSummary', ...
        'Run run_baseline_comparison.m first. Missing file: %s', summaryFile);
end

summaryData = load(summaryFile);
summaryRows = summaryData.summaryRows;
config = build_config(summaryData.configSummary);
targetOrbit = config.targetOrbit;

records = [load_proposed_records(resultsDir), load_result_records(summaryRows)];
if isempty(records)
    error('draw_baseline_comparison:NoResults', ...
        'No successful baseline result files found in %s.', resultsDir);
end

[groupKeys, groupLabels, groupIndex] = unique_groups(records);
colors = lines(numel(groupKeys));

%% Figure 1: Best-so-far loss vs evaluations
figure
hold on
grid on
box on

for groupIdx = 1:numel(groupKeys)
    groupRecords = records(groupIndex == groupIdx);
    curves = collect_curves(groupRecords);
    evalAxis = 1:size(curves, 2);
    medianCurve = row_nan_percentile(curves, 50);
    q1Curve = row_nan_percentile(curves, 25);
    q3Curve = row_nan_percentile(curves, 75);

    semilogy(evalAxis, medianCurve, 'LineWidth', 1.5, ...
        'Color', colors(groupIdx, :), 'DisplayName', groupLabels{groupIdx});
    semilogy(evalAxis, q1Curve, ':', 'LineWidth', 0.75, ...
        'Color', colors(groupIdx, :), 'HandleVisibility', 'off');
    semilogy(evalAxis, q3Curve, ':', 'LineWidth', 0.75, ...
        'Color', colors(groupIdx, :), 'HandleVisibility', 'off');
end

xlabel('Forward evaluations')
ylabel('Best-so-far loss')
title('Best-so-far loss vs evaluations')
legend('Location', 'best')
set(gca, 'YScale', 'log')

%% Figure 2: Final best loss across seeds
figure
hold on
grid on
box on

for groupIdx = 1:numel(groupKeys)
    groupRecords = records(groupIndex == groupIdx);
    finalLosses = [groupRecords.finalLoss].';
    draw_distribution_box(groupIdx, finalLosses, colors(groupIdx, :));
end

xlim([0.5, numel(groupKeys) + 0.5])
set(gca, 'XTick', 1:numel(groupKeys), 'XTickLabel', groupLabels, ...
    'XTickLabelRotation', 30, 'YScale', 'log')
ylabel('Final best loss')
title('Final best loss across seeds')

%% Figure 3: Runtime vs final best loss
figure
hold on
grid on
box on

for groupIdx = 1:numel(groupKeys)
    groupRecords = records(groupIndex == groupIdx);
    runtimes = [groupRecords.runtime].';
    finalLosses = [groupRecords.finalLoss].';
    scatter(runtimes, finalLosses, 36, colors(groupIdx, :), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'DisplayName', groupLabels{groupIdx});
end

xlabel('Runtime (s)')
ylabel('Final best loss')
title('Runtime vs final best loss')
set(gca, 'YScale', 'log')
legend('Location', 'best')

%% Figure 4: Best parameter time series comparison
figure
stateCount = size(targetOrbit.y, 2);
representatives = select_group_representatives(records, groupIndex, numel(groupKeys));
plotPhase = linspace(0, 1, config.lossOptions.compareNumPoints + 1).';
plotPhase(end) = [];
targetY = resample_orbit_y(targetOrbit, plotPhase);

for stateIdx = 1:stateCount
    subplot(stateCount, 1, stateIdx)
    hold on
    grid on
    box on

    plot(plotPhase, targetY(:, stateIdx), 'k-', 'LineWidth', 1.8, ...
        'DisplayName', 'target');

    for groupIdx = 1:numel(groupKeys)
        record = representatives(groupIdx);
        if isempty(record.bestOrbit)
            continue
        end

        candidateY = resample_orbit_y(record.bestOrbit, plotPhase);
        candidateY = align_candidate_phase(candidateY, targetY);
        plot(plotPhase, candidateY(:, stateIdx), 'LineWidth', 1.1, ...
            'Color', colors(groupIdx, :), 'DisplayName', groupLabels{groupIdx});
    end

    ylabel(sprintf('State %d', stateIdx))
    if stateIdx == 1
        title('Best parameter time series comparison')
        legend('Location', 'bestoutside')
    end
    if stateIdx == stateCount
        xlabel('Normalized phase')
    end
end

%% Local functions
function records = load_result_records(summaryRows)
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
    if ~strcmp(row.status, 'success') || isempty(row.resultFile) || ~isfile(row.resultFile)
        continue
    end

    fileData = load(row.resultFile);
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

function records = load_proposed_records(resultsDir)
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

proposedFile = fullfile(resultsDir, 'proposed_method', 'proposed_method_result.mat');
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

function draw_distribution_box(x, values, color)
values = values(isfinite(values));
if isempty(values)
    return
end

q1 = percentile_value(values, 25);
q2 = percentile_value(values, 50);
q3 = percentile_value(values, 75);
vMin = min(values);
vMax = max(values);
boxWidth = 0.28;

patch(x + boxWidth * [-1 1 1 -1], [q1 q1 q3 q3], color, ...
    'FaceAlpha', 0.18, 'EdgeColor', color, 'LineWidth', 1.0);
plot([x - boxWidth, x + boxWidth], [q2 q2], '-', 'Color', color, ...
    'LineWidth', 1.5);
plot([x x], [vMin q1], '-', 'Color', color, 'LineWidth', 0.8);
plot([x x], [q3 vMax], '-', 'Color', color, 'LineWidth', 0.8);
plot([x - boxWidth / 2, x + boxWidth / 2], [vMin vMin], '-', ...
    'Color', color, 'LineWidth', 0.8);
plot([x - boxWidth / 2, x + boxWidth / 2], [vMax vMax], '-', ...
    'Color', color, 'LineWidth', 0.8);

jitter = linspace(-0.08, 0.08, numel(values)).';
if isscalar(values)
    jitter = 0;
end
scatter(x + jitter, values, 24, color, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
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
