% Plot the FHN FMAM amplitude-invariant sweep against the normal-form FM curve.

ensure_paths();

scriptDir = fileparts(mfilename('fullpath'));
resultFile = fullfile(scriptDir, 'fhn_fmam_results.mat');

loaded = load(resultFile, 'result');
result = loaded.result;
normalData = build_fhn_normal_form_data(struct('assert_match', false));

fmamParams = vertcat(result.fm.params);
f12Span = linspace(-1, 0.4, 200);

fig = figure('Name', 'FHN FMAM amplitude-invariant sweep');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
hold(ax, 'on');
plot(ax, normalData.curves.fm.f11, normalData.curves.fm.f12, ...
    '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.6);
plot(ax, fmamParams(:, 5), fmamParams(:, 6), ...
    's', 'Color', [0.6500, 0.1800, 0.0500], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.1);
grid(ax, 'on');
xlabel(ax, 'f_{11}');
ylabel(ax, 'f_{12}');
title(ax, 'Amplitude invariant');
legend(ax, {'normal-form FM', 'FMAM FM'}, 'Location', 'best');
xlim(ax, [-0.2, 0.4]);
ylim(ax, [-1.5, 1]);

ax = nexttile(layout);
hold(ax, 'on');
plot(ax, f12Span, normalData.core.a21_over_a12 .* f12Span, ...
    '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.4);
plot(ax, fmamParams(:, 6), fmamParams(:, 7), ...
    's', 'Color', [0.6500, 0.1800, 0.0500], ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.1);
grid(ax, 'on');
xlabel(ax, 'f_{12}');
ylabel(ax, 'f_{21}');
title(ax, 'Structure constraint');
legend(ax, {'scale-invariance condition', 'FMAM FM'}, 'Location', 'best');
xlim(ax, [-1, 0.4]);
ylim(ax, [-0.05, 0.15]);
