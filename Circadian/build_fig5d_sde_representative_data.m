clear
clc

%% Paths
disp('%% Paths');
scriptDir = fileparts(mfilename('fullpath'));
markerDataDir = fullfile(scriptDir, 'data', 'fig5d', 'markers');
sdeDataDir = fullfile(scriptDir, 'data', 'fig5d', 'sde');
representativeDataFile = fullfile(sdeDataDir, 'fig5d_sde_representative.mat');

addpath(fullfile(scriptDir, 'SDE_simulation'));
mkdir(sdeDataDir);

%% Files
disp('%% Files');
markerMaximum = 0.12;
markerPeriods = [23.5, 24.0, 24.5];

markerFiles = cell(1, numel(markerPeriods));
for i = 1:numel(markerPeriods)
    markerFiles{i} = fullfile(markerDataDir, ...
        sprintf('fig5d_marker_M%s_T%s.mat', maximum_tag(markerMaximum), period_tag(markerPeriods(i))));
end

%% Simulation settings
disp('%% Simulation settings');
representativeSeeds = [0, 0, 0];
warmupT = 1000;
measureT = 72;

representativeOptions = struct();
representativeOptions.warmupT = warmupT;
representativeOptions.measureT = measureT;
representativeOptions.T = warmupT + measureT;
representativeOptions.dt = 0.005;
representativeOptions.sigma = 0.1;
representativeOptions.noiseClass = {'o', 'o', 'o'};

%% Simulate representative traces
disp('%% Simulate representative traces');
representative = repmat(struct( ...
    'targetPeriod', [], ...
    'seed', [], ...
    'Parameters', [], ...
    't', [], ...
    'X', []), 1, numel(markerFiles));

for i = 1:numel(markerFiles)
    fprintf('Representative marker %d/%d: T = %.1f, seed = %d\n', ...
        i, numel(markerFiles), markerPeriods(i), representativeSeeds(i));
    loaded = load(markerFiles{i}, 'Parameters', 'TS', 'period');

    y = loaded.TS{2};
    y0 = y(1, :).';

    runOptions = representativeOptions;
    runOptions.seed = representativeSeeds(i);

    sim = Circadian_SDE(loaded.Parameters, y0, runOptions);
    measureMask = sim.t >= representativeOptions.warmupT;

    representative(i).targetPeriod = loaded.period(1);
    representative(i).seed = representativeSeeds(i);
    representative(i).Parameters = reshape(loaded.Parameters, 1, []);
    representative(i).t = sim.t(measureMask) - representativeOptions.warmupT;
    representative(i).X = sim.X(measureMask, :);
end

%% Save data
disp('%% Save data');
save(representativeDataFile, ...
    'representativeSeeds', ...
    'representativeOptions', ...
    'representative');

fprintf('Saved data: %s\n', representativeDataFile);

function tag = maximum_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end
