clear
clc
%%
netName = 'BA';
dynamicName = 'GRN';
propname = 'amp';
weight_per = 0.03;
N = 100;
sourceSequenceIndex = 2;
threshold = 0.05;
n_target = 50;
scatterExportResolution = 600;

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
add_cbrewer_path(repoDir);
dataFile = fullfile( ...
    scriptDir, 'Ergodic data', sprintf('N = %d', N), ...
    sprintf('%s (net = %s, source_seq = %d, weight_per = %s).mat', ...
    dynamicName, netName, sourceSequenceIndex, num2str(weight_per)));
figDir = fullfile(scriptDir, 'temp_fig');
ensure_directory(figDir);

data = networkexp.load_ergodic_results(dataFile);
prop_ori = get_primary_observable(data.baseline.observables, propname, N);
sourceNodeOrder = data.sourceNodeSequence.nodeOrder;

proportion_mat = NaN(N, N);
for nTarget = 1:N
    targetNodes = sourceNodeOrder(1:nTarget);
    for idxRepeat = 1:numel(data.repeats)
        samples = data.repeats(idxRepeat).plottableSamples;
        if isempty(samples)
            continue
        end
        proportion = zeros(1, numel(samples));
        for i = 1:numel(samples)
            sample_i = get_primary_observable(samples(i).observables, propname, N);
            proportion(i) = sum( ...
                ((sample_i(targetNodes) - prop_ori(targetNodes)) ./ prop_ori(targetNodes)) > threshold) / nTarget;
        end
        proportion_mat(idxRepeat, nTarget) = mean(proportion);
    end
end

fig = figure;
imagesc(1:N, 1:N, proportion_mat);
set(gca, 'YDir', 'reverse')
hold on
plot([1, N], [1, N], '--', 'Color', 'w', 'LineWidth', 1)
hold off
colormap(gca, get_heatmap_colormap(256))
clim([0, 1])
colorbarHandle = colorbar;
xlabel('Number of targets');
ylabel('Number of edges modulated');
set(gca, ...
    'XTick', [1, N], ...
    'YTick', [1, N], ...
    'XTickLabel', {'1', num2str(N)}, ...
    'YTickLabel', {'1', num2str(N)}, ...
    'TickLength', [0, 0], ...
    'fontsize', 8)
ax = gca;
ax.Position = [0.18, 0.20, 0.62, 0.72];
colorbarHandle.Position = [0.83, 0.20, 0.035, 0.72];
ax.XLabel.Units = 'normalized';
ax.YLabel.Units = 'normalized';
ax.XLabel.Position(2) = -0.055;
ax.YLabel.Position(1) = -0.055;
box on
add_axis_arrows(gca)
fig.PaperUnits = "centimeters";
width = 6;
height = 6 * 0.75;
fig.PaperSize = [width, height];
fig.PaperPosition = [0, 0, width, height];
figname = ['Heatmap_' dynamicName '_source_seq = ' num2str(sourceSequenceIndex) '_' netName];
print(fig, fullfile(figDir, [figname '.pdf']), '-dpdf');

targetNodes = sourceNodeOrder(1:n_target);
numScatterPoints = 0;
for idxRepeat = 1:numel(data.repeats)
    numScatterPoints = numScatterPoints + numel(data.repeats(idxRepeat).plottableSamples);
end
proportion_sample = NaN(numScatterPoints, 2);
sampleCursor = 0;
for idxRepeat = 1:numel(data.repeats)
    samples = data.repeats(idxRepeat).plottableSamples;
    for i = 1:numel(samples)
        sample_i = get_primary_observable(samples(i).observables, propname, N);
        modulationRate = sum( ...
            ((sample_i(targetNodes) - prop_ori(targetNodes)) ./ prop_ori(targetNodes)) > threshold) / n_target;
        sampleCursor = sampleCursor + 1;
        proportion_sample(sampleCursor, :) = [idxRepeat, modulationRate];
    end
end
proportion_sample = proportion_sample(1:sampleCursor, :);

figure
if isempty(proportion_sample)
    error('draw_ergodic:NoPlottableSamples', ...
        'No successful or candidate periodic samples were found in %s.', dataFile);
end
scatter(proportion_sample(:, 1), proportion_sample(:, 2), 3, 'filled', 'MarkerFaceAlpha', 0.5);
title(['Scatter of ' netName ' ' dynamicName])
xlabel('Number of edges modulated')
ylabel('Modulation rate')
axis([1 N 0 1])
xline(n_target, 'r--', 'LineWidth', 1)
box on
set(gca, 'fontsize', 8)

fig = gcf;
fig.PaperUnits = "centimeters";
width = 6;
height = 6 * 0.75;
fig.PaperSize = [width, height];
fig.PaperPosition = [0, 0, width, height];
figname = ['Scatter_' dynamicName '_n_target = ' num2str(n_target) ...
    '_source_seq = ' num2str(sourceSequenceIndex) '_' netName];
save_scatter_figure(fig, fullfile(figDir, [figname '.pdf']), scatterExportResolution);

function values = get_primary_observable(observables, propname, nodeCount)
switch lower(propname)
    case {'amp', 'amplitude'}
        values = observables.amplitude(1:nodeCount);
    case {'max', 'maxvariable'}
        values = observables.maxVariable(1:nodeCount);
    case {'min', 'minvariable'}
        values = observables.minVariable(1:nodeCount);
    otherwise
        error('draw_ergodic:UnsupportedObservable', ...
            'Observable "%s" is not supported for target-wise modulation plots.', propname);
end
end

function ensure_directory(folderPath)
if exist(folderPath, 'dir') ~= 7
    mkdir(folderPath);
end
end

function add_cbrewer_path(repoDir)
cbrewerDir = fullfile(repoDir, '+utils', 'cbrewer2');
if exist(cbrewerDir, 'dir') == 7 && ~contains(path, cbrewerDir)
    addpath(cbrewerDir);
end
end

function add_axis_arrows(ax)
fig = ancestor(ax, 'figure');
drawnow
axPosition = ax.Position;
annotation(fig, 'arrow', ...
    [axPosition(1), axPosition(1) + axPosition(3)], ...
    [axPosition(2) - 0.085, axPosition(2) - 0.085], ...
    'LineWidth', 0.75, 'Color', 'k');
annotation(fig, 'arrow', ...
    [axPosition(1) - 0.085, axPosition(1) - 0.085], ...
    [axPosition(2) + axPosition(4), axPosition(2)], ...
    'LineWidth', 0.75, 'Color', 'k');
end

function save_scatter_figure(fig, filePath, resolution)
exportgraphics(fig, filePath, 'ContentType', 'image', 'Resolution', resolution);
end

function cmap = get_heatmap_colormap(numColors)
if nargin < 1
    numColors = 256;
end
if exist('cbrewer2', 'file') ~= 2
    error('draw_ergodic:MissingCbrewer2', ...
        'Could not find cbrewer2. Expected +utils/cbrewer2/cbrewer2.m under the repository root.');
end
cmap = flipud(cbrewer2('RdYlBu', numColors, 'linear', 'rgb'));
end
