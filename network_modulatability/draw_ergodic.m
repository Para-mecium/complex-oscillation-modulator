clear
clc

netName = 'ER';
dynamicName = 'FHN';
propname = 'amp';
weight_per = 0.3;
threshold = 0.05;
n_target = 50;

scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile( ...
    scriptDir, 'Ergodic data', ...
    [dynamicName ' (net = ' netName ', weight_per = ' num2str(weight_per) ')' '.mat']);
figDir = fullfile(scriptDir, 'temp_fig');
ensure_directory(figDir);

data = networkexp.load_ergodic_results(dataFile);
N = data.nodeCount;
prop_ori = get_primary_observable(data.baseline.observables, propname, N);

proportion_mat = NaN(N, N);
for nTarget = 1:N
    for idxRepeat = 1:numel(data.repeats)
        samples = data.repeats(idxRepeat).plottableSamples;
        if isempty(samples)
            continue
        end
        proportion = zeros(1, numel(samples));
        for i = 1:numel(samples)
            sample_i = get_primary_observable(samples(i).observables, propname, N);
            proportion(i) = sum( ...
                ((sample_i(1:nTarget) - prop_ori(1:nTarget)) ./ prop_ori(1:nTarget)) > threshold) / nTarget;
        end
        proportion_mat(idxRepeat, nTarget) = mean(proportion);
    end
end

figure
h = heatmap(proportion_mat);
h.GridVisible = false;
h.XDisplayData = 1:N;
h.YDisplayData = 1:N;
title(['Heatmap of ' netName ' ' dynamicName])
ylabel('Number of edges modulated');
xlabel('Number of targets');

proportion_sample = [];
for idxRepeat = 1:numel(data.repeats)
    samples = data.repeats(idxRepeat).plottableSamples;
    for i = 1:numel(samples)
        sample_i = get_primary_observable(samples(i).observables, propname, N);
        modulationRate = sum( ...
            ((sample_i(1:n_target) - prop_ori(1:n_target)) ./ prop_ori(1:n_target)) > threshold) / n_target;
        proportion_sample = [proportion_sample; idxRepeat, modulationRate]; %#ok<AGROW>
    end
end

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
box on
set(gca, 'fontsize', 8)

fig = gcf;
fig.PaperUnits = "centimeters";
width = 6;
height = 6 * 0.75;
fig.PaperSize = [width, height];
fig.PaperPosition = [0, 0, width, height];
figname = ['Scatter_' dynamicName '_n_target = ' num2str(n_target) '_' netName];
print(fig, fullfile(figDir, [figname '.pdf']), '-dpdf');

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
