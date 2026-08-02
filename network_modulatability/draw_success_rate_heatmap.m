clear
clc
%%
netName = 'ER';
dynamicName = 'FHN';
N = 100;
sourceSequenceIndex = 1;
weight_per = 0.3;

scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile( ...
    scriptDir, 'Ergodic data', sprintf('N = %d', N), ...
    sprintf('%s (net = %s, source_seq = %d, weight_per = %s).mat', ...
    dynamicName, netName, sourceSequenceIndex, num2str(weight_per)));
figDir = fullfile(scriptDir, 'temp_fig');
ensure_directory(figDir);

data = networkexp.load_ergodic_results(dataFile);

successRateVec = NaN(1, N);
for idxRepeat = 1:numel(data.repeats)
    rate = (data.repeats(idxRepeat).successCount + data.repeats(idxRepeat).candidateCount) ...
        / data.repeats(idxRepeat).attemptCount;
    successRateVec(idxRepeat) = 1-rate;
end

fig = figure;
plot(1:N, successRateVec, 'k-', 'LineWidth', 1)
xlabel('Number of edges modulated');
ylabel('Exclusion rate');
set(gca, ...
    'XTick', 0:N/5:N, ...
    'fontsize', 8)
axis([0 N 0 1])
ax = gca;
ax.Position = [0.18, 0.20, 0.72, 0.72];
ax.XLabel.Units = 'normalized';
ax.YLabel.Units = 'normalized';
ax.XLabel.Position(2) = -0.055;
ax.YLabel.Position(1) = -0.055;
grid on
box on
fig.PaperUnits = "centimeters";
width = 6;
height = 6 * 0.75;
fig.PaperSize = [width, height];
fig.PaperPosition = [0, 0, width, height];
figname = ['SuccessRate_' dynamicName '_source_seq = ' num2str(sourceSequenceIndex) '_' netName];
print(fig, fullfile(figDir, [figname '.pdf']), '-dpdf');

function ensure_directory(folderPath)
if exist(folderPath, 'dir') ~= 7
    mkdir(folderPath);
end
end
