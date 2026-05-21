% Plot the cancer FMAM period-invariant sweep against the normal-form AM curve.

scriptDir = fileparts(mfilename('fullpath'));
normalFormDir = fileparts(scriptDir);
repoDir = fileparts(normalFormDir);
addpath(repoDir);
addpath(normalFormDir);
addpath(scriptDir);

resultFile = fullfile(scriptDir, 'cancer_fmam_results.mat');

loaded = load(resultFile, 'result');
result = loaded.result;
cancerOpts = struct( ...
    'assert_match', false, ...
    'epsilon', result.baseline.params(5));
normalData = build_cancer_normal_form_data(cancerOpts);

fmamParams = vertcat(result.am.params);
f12Span = linspace(-4, 1, 200);

fig = figure('Name', 'Cancer FMAM period-invariant sweep');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
hold(ax, 'on');
plot(ax, normalData.curves.am.f11, normalData.curves.am.f12, ...
    '-', 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1.6);
plot(ax, fmamParams(:, 6), fmamParams(:, 7), ...
    'o', 'Color', [0.0000, 0.2500, 0.6000], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.1);
grid(ax, 'on');
xlabel(ax, 'f_{11}');
ylabel(ax, 'f_{12}');
title(ax, 'Period invariant');
legend(ax, {'normal-form AM', 'FMAM AM'}, 'Location', 'best');
xlim(ax, [-2, 3]);
ylim(ax, [-4, 1]);

ax = nexttile(layout);
hold(ax, 'on');
plot(ax, f12Span, normalData.core.a21_over_a12 .* f12Span, ...
    '-', 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1.4);
plot(ax, fmamParams(:, 7), fmamParams(:, 8), ...
    'o', 'Color', [0.0000, 0.2500, 0.6000], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.1);
grid(ax, 'on');
xlabel(ax, 'f_{12}');
ylabel(ax, 'f_{21}');
title(ax, 'Structure constraint');
legend(ax, {'scale-invariance condition', 'FMAM AM'}, 'Location', 'best');
xlim(ax, [-4, 1]);
ylim(ax, [-0.5, 2]);
