clear
clc

scriptDir = fileparts(mfilename('fullpath'));
data = load(fullfile(scriptDir, 'figS17_repressilator_data.mat'));

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4, 4, 13, 4.6]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:numel(data.phaseResults)
    result = data.phaseResults(i);
    t = result.TS{1};
    y = result.TS{2};
    [~, idx1] = max(y(:, 4));
    t = t - t(idx1);
    t(t < 0) = t(t < 0) + result.period;
    [t, order] = sort(t);
    y = y(order, :);

    nexttile
    hold on
    plot(t, y(:, 4), 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1.5)
    plot(t, y(:, 5), 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.5)
    plot(result.phase12, max(y(:, 5)), 'k.', 'MarkerSize', 10)
    text(result.phase12 + 1.2, max(y(:, 5)), sprintf('t* = %.0f', result.targetPhase), ...
        'FontSize', 8, 'VerticalAlignment', 'middle')
    xlim([0, 50])
    ylim([0, 150])
    grid on
    box on
    xlabel('Time (min)', 'FontSize', 9)
    ylabel('Concentration (a.u.)', 'FontSize', 9)
    set(gca, 'FontSize', 8, 'LineWidth', 0.8)
    text(0.02, 1.07, char('A' + i - 1), 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 11)
end

exportgraphics(fig, fullfile(scriptDir, 'Fig_S17.png'), 'Resolution', 300)
exportgraphics(fig, fullfile(scriptDir, 'Fig_S17.pdf'), 'ContentType', 'vector')
