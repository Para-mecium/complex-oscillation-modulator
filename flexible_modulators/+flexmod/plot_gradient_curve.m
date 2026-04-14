function h = plot_gradient_curve(ax, x, y, colorValue, lineWidth)
if nargin < 5
    lineWidth = 2.5;
end

x = x(:).';
y = y(:).';
colorValue = colorValue(:).';

if numel(x) < 2
    h = plot(ax, x, y, 'LineWidth', lineWidth);
    return
end

h = surface(ax, [x; x], [y; y], zeros(2, numel(x)), [colorValue; colorValue], ...
    'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', lineWidth);
end
