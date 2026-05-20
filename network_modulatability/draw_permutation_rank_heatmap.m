clear
clc

%% User settings
N = 100;
dynamicName = 'GRN';
netName = 'BA';
sourceSequenceIndex = 4;
weight_per = 0.03;

% Options: 'sourceSequence' or 'identity'
rankMode = 'sourceSequence';

% Set to [] to skip target highlighting.
nTarget = [];

scriptDir = fileparts(mfilename('fullpath'));
figDir = fullfile(scriptDir, 'temp_fig');
ensure_directory(figDir);

[nodeOrder, labelText] = load_rank_source( ...
    scriptDir, dynamicName, netName, N, sourceSequenceIndex, weight_per, rankMode);
nodeRank = order_to_rank(nodeOrder, N);

fig = figure;
imagesc(1:N, 1, nodeRank);
set(gca, 'YDir', 'normal')
colormap(gca, flipud(gray_scale(256, 0.1, 0.9)))
clim([1, N])
colorbarHandle = colorbar;
colorbarHandle.Label.String = 'Rank in permutation';
hold on
if ~isempty(nTarget)
    targetNodes = nodeOrder(1:min(nTarget, N));
    scatter(targetNodes, ones(size(targetNodes)), 9, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'w', ...
        'MarkerEdgeAlpha', 0.9, 'LineWidth', 0.35);
end
hold off

xlabel('Node id');
yticks([])
% yticklabels({labelText})
set(gca, ...
    'XTick', get_axis_ticks(N), ...
    'TickLength', [0, 0], ...
    'fontsize', 8)
ax = gca;
ax.Position = [0.12, 0.34, 0.70, 0.28];
colorbarHandle.Position = [0.85, 0.34, 0.035, 0.28];
box on

fig.PaperUnits = "centimeters";
width = 9;
height = 2.4;
fig.PaperSize = [width, height];
fig.PaperPosition = [0, 0, width, height];

fileToken = sprintf('%s_%s_N%d_%s', dynamicName, netName, N, rankMode);
if strcmpi(rankMode, 'sourceSequence')
    fileToken = sprintf('%s_source_seq_%d', fileToken, sourceSequenceIndex);
end
if ~isempty(nTarget)
    fileToken = sprintf('%s_n_target_%d', fileToken, nTarget);
end
print(fig, fullfile(figDir, ['Permutation_rank_' fileToken '.pdf']), '-dpdf');

function [nodeOrder, labelText] = load_rank_source( ...
        scriptDir, dynamicName, netName, N, sourceSequenceIndex, weightPer, rankMode)
switch lower(char(string(rankMode)))
    case 'identity'
        nodeOrder = 1:N;
        labelText = 'Original rank';

    case 'sourcesequence'
        dataFile = fullfile( ...
            scriptDir, 'Ergodic data', sprintf('N = %d', N), ...
            sprintf('%s (net = %s, source_seq = %d, weight_per = %s).mat', ...
            dynamicName, netName, sourceSequenceIndex, num2str(weightPer)));
        data = load(dataFile, 'sourceNodeSequence');
        nodeOrder = reshape(data.sourceNodeSequence.nodeOrder, 1, []);
        labelText = sprintf('Source seq %d', sourceSequenceIndex);
end
end

function nodeRank = order_to_rank(nodeOrder, N)
nodeRank = zeros(1, N);
nodeRank(nodeOrder) = 1:N;
end

function ticks = get_axis_ticks(N)
if N <= 25
    ticks = 1:N;
else
    ticks = unique([1, round(N / 4), round(N / 2), round(3 * N / 4), N]);
end
end

function ensure_directory(folderPath)
if exist(folderPath, 'dir') ~= 7
    mkdir(folderPath);
end
end

function cmap = gray_scale(numColors, lowValue, highValue)
values = linspace(lowValue, highValue, numColors).';
cmap = [values, values, values];
end
