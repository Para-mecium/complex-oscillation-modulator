clear
clc
%%
dynamicName = 'FHN';
netName = 'ER';
N = 100;
n_per = 50;
mksize = 0.5;
arsize = 2;
lnwidth = 0.2;
modulatedEdgeAlpha = 0.5;

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

G = build_plot_graph(baseline.weightMatrix);
figure
net = plot(G, "Layout", "force", NodeLabel={}, MarkerSize=mksize, ArrowSize=arsize, LineWidth=lnwidth);
net.XData = net.XData * 2;
net.YData = net.YData * 2;
axis off
box off
nodePositions = [net.XData', net.YData'];

format_figure(gcf)
save_network_figure(gcf, fullfile(figDir, [netName '_origin.pdf']));

G = build_plot_graph(baseline.weightMatrix);
figure
net_per = plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), ...
    NodeLabel={}, MarkerSize=mksize, ArrowSize=arsize, LineWidth=lnwidth);
hold on
plot_perturbed_edges(perturbed, N, nodePositions, arsize, lnwidth * 3, modulatedEdgeAlpha);
hold off
axis off
box off

format_figure(gcf)
save_network_figure(gcf, fullfile(figDir, [netName '_n_per = ' num2str(n_per) '.pdf']));

function G = build_plot_graph(weightMatrix)
adjacency = spones(sparse(weightMatrix));
G = digraph(adjacency);
end

function save_network_figure(fig, filePath)
print(fig, filePath, '-dpdf', '-painters');
end

function plot_perturbed_edges(perturbed, nodeCount, nodePositions, arrowSize, lineWidth, edgeAlpha)
edgeList = perturbed.perturbationEdges;
if isempty(edgeList)
    return
end

edgeList = edgeList(all(edgeList >= 1 & edgeList <= nodeCount, 2), :);
if isempty(edgeList)
    return
end

edgeGraph = digraph(edgeList(:, 1), edgeList(:, 2), [], nodeCount);
edgePlot = plot(edgeGraph, ...
    'XData', nodePositions(:, 1), ...
    'YData', nodePositions(:, 2), ...
    'NodeLabel', {}, ...
    'Marker', 'none', ...
    'EdgeColor', 'r', ...
    'EdgeAlpha', edgeAlpha, ...
    'ArrowSize', arrowSize, ...
    'LineWidth', lineWidth);
uistack(edgePlot, 'top');
end

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
