clear
clc

%% Paths
disp('%% Paths');
scriptDir = fileparts(mfilename('fullpath'));
markerDataDir = fullfile(scriptDir, 'data', 'fig5d', 'markers');
sdeDataDir = fullfile(scriptDir, 'data', 'fig5d', 'sde');
repeatDataFile = fullfile(sdeDataDir, 'fig5d_sde_repeat.mat');

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
repeatSeeds = 1:100;
warmupT = 1000;
psdMeasureT = 10000;
distributionMeasureT = 24;
psdBandPercentiles = [10, 90];

repeatOptions = struct();
repeatOptions.T = warmupT + psdMeasureT;
repeatOptions.dt = 0.01;
repeatOptions.sigma = 0.1;
repeatOptions.noiseClass = {'o', 'o', 'o'};
repeatOptions.psdObservable = 'Ptot';
psdMeasureCount = round(psdMeasureT / repeatOptions.dt) + 1;
distributionMeasureCount = round(distributionMeasureT / repeatOptions.dt) + 1;

%% Simulate repeat batches
disp('%% Simulate repeat batches');
psdStats = repmat(struct( ...
    'targetPeriod', [], ...
    'Parameters', [], ...
    'frequency', [], ...
    'mean', [], ...
    'lower', [], ...
    'upper', []), 1, numel(markerFiles));
distributionStats = repmat(struct( ...
    'targetPeriod', [], ...
    'Parameters', [], ...
    'maxPtot', []), 1, numel(markerFiles));

for i = 1:numel(markerFiles)
    fprintf('Repeat marker %d/%d: T = %.1f\n', ...
        i, numel(markerFiles), markerPeriods(i));
    loaded = load(markerFiles{i}, 'Parameters', 'TS', 'period');

    y = loaded.TS{2};
    y0 = y(1, :).';
    parameters = loaded.Parameters;
    targetPeriod = loaded.period(1);

    psdValues = zeros(numel(repeatSeeds), floor(psdMeasureCount / 2) + 1);
    maxPtot = zeros(numel(repeatSeeds), 1);
    frequencyVectors = cell(1, numel(repeatSeeds));

    parfor j = 1:numel(repeatSeeds)
        fprintf('  Repeat %d/%d: seed = %d\n', ...
            j, numel(repeatSeeds), repeatSeeds(j));
        runOptions = repeatOptions;
        runOptions.seed = repeatSeeds(j);

        sim = Circadian_SDE(parameters, y0, runOptions);
        measureMask = sim.t >= warmupT;
        X = sim.X(measureMask, :);

        psdX = X(1:psdMeasureCount, :);
        Ptot = psdX(:, 2) + psdX(:, 3);
        [frequency, spectrum] = single_sided_periodogram(Ptot, repeatOptions.dt);
        frequencyVectors{j} = frequency(:);
        psdValues(j, :) = spectrum(:).';

        distributionX = X(1:distributionMeasureCount, :);
        Ptot = distributionX(:, 2) + distributionX(:, 3);
        maxPtot(j) = max(Ptot);
    end

    psdBand = prctile(psdValues, psdBandPercentiles, 1);

    psdStats(i).targetPeriod = targetPeriod;
    psdStats(i).Parameters = reshape(parameters, 1, []);
    psdStats(i).frequency = frequencyVectors{1};
    psdStats(i).mean = mean(psdValues, 1).';
    psdStats(i).lower = psdBand(1, :).';
    psdStats(i).upper = psdBand(2, :).';

    distributionStats(i).targetPeriod = targetPeriod;
    distributionStats(i).Parameters = reshape(parameters, 1, []);
    distributionStats(i).maxPtot = maxPtot;
end

%% Save data
disp('%% Save data');
save(repeatDataFile, ...
    'repeatSeeds', ...
    'repeatOptions', ...
    'warmupT', ...
    'psdMeasureT', ...
    'distributionMeasureT', ...
    'psdBandPercentiles', ...
    'psdStats', ...
    'distributionStats');

fprintf('Saved data: %s\n', repeatDataFile);

function [frequency, spectrum] = single_sided_periodogram(samples, dt)
samples = reshape(double(samples), [], 1);
samples = samples - mean(samples);

n = numel(samples);
fs = 1 / dt;
fftValues = fft(samples);
powerValues = abs(fftValues).^2 / (fs * max(n, 1));

halfIdx = floor(n / 2) + 1;
frequency = (0:halfIdx - 1).' * (fs / n);
spectrum = powerValues(1:halfIdx);
if numel(spectrum) > 2
    spectrum(2:end-1) = 2 * spectrum(2:end-1);
end
end

function tag = maximum_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end
