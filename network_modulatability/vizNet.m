dynamicName = 'FHN';
netName = 'ER';
N = 100;
n_per = 50;
mksize = 0.5;
arsize = 2;
lnwidth = 0.2;

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

G = digraph(baseline.weightMatrix);
figure
net = plot(G, "Layout", "force", NodeLabel={}, MarkerSize=mksize, ArrowSize=arsize, LineWidth=lnwidth);
net.XData = net.XData * 2;
net.YData = net.YData * 2;
nodePositions = [net.XData', net.YData'];

format_figure(gcf)
print(gcf, fullfile(figDir, [netName '_origin.pdf']), '-dpdf');

G = digraph(perturbed.weightMatrix);
figure
net_per = plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), ...
    NodeLabel={}, MarkerSize=mksize, ArrowSize=arsize, LineWidth=lnwidth);
if ~isempty(perturbed.edgeIndices)
    highlight(net_per, 'Edges', perturbed.edgeIndices, 'EdgeColor', 'r', 'LineWidth', lnwidth * 3)
end

format_figure(gcf)
print(gcf, fullfile(figDir, [netName '_n_per = ' num2str(n_per) '.pdf']), '-dpdf');

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
