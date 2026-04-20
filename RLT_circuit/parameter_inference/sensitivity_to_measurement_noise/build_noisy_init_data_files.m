clear
clc

noiseLevels = [0.01 0.02 0.05];
baseSeed = 1;
N = 100;

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
outputDir = fullfile(scriptDir, 'noisy_init_data_files');

sourceData = load(fullfile(parameterInferenceDir, 'initData_circuit.mat'));
t = sourceData.t(:);
y = sourceData.y;

periodBase = t(end) - t(1);
tTwoPeriod = [t; t(2:end) + periodBase];
yTwoPeriod = [y; y(2:end, :)];

totalDuration = tTwoPeriod(end) - tTwoPeriod(1);
medianDt = median(diff(tTwoPeriod));
frequencyGrid = linspace(1 / totalDuration, 1 / (2 * medianDt), 4000).';

if ~isfolder(outputDir)
    mkdir(outputDir);
end

for noiseIdx = 1:numel(noiseLevels)
    noiseLevel = noiseLevels(noiseIdx);
    noiseLevelDirName = ['noise_level_' strrep(num2str(noiseLevel), '.', 'p')];
    noiseLevelDir = fullfile(outputDir, noiseLevelDirName);

    if ~isfolder(noiseLevelDir)
        mkdir(noiseLevelDir);
    end

    for i = 1:N
        seed = baseSeed + (i - 1);
        rng(seed, 'twister');

        yNoisy = yTwoPeriod + noiseLevel * randn(size(yTwoPeriod));
        TS = {tTwoPeriod, yNoisy};

        varMax = reshape(max(yNoisy, [], 1), 1, []);
        varMin = reshape(min(yNoisy, [], 1), 1, []);
        varAmp = reshape((varMax - varMin)/2, 1, []);
        period = estimate_period_lomb(tTwoPeriod, yNoisy, frequencyGrid);

        saveFile = fullfile(noiseLevelDir, sprintf('initData_%03d.mat', i));
        save(saveFile, 'TS', 'seed', 'noiseLevel', 'period', 'varAmp', 'varMax', 'varMin');
        msg = sprintf('Saved noiseLevel=%.6g, seed=%d, period=%.6g: %s', ...
            noiseLevel, seed, period, saveFile);
        fprintf('%s\n', msg);
    end
end

function period = estimate_period_lomb(t, y, frequencyGrid)
peakPower = zeros(1, size(y, 2));
peakFrequency = zeros(1, size(y, 2));

for stateIdx = 1:size(y, 2)
    centeredSignal = y(:, stateIdx) - mean(y(:, stateIdx));
    [powerSpectrum, frequency] = plomb(centeredSignal, t, frequencyGrid, 'psd');
    [peakPower(stateIdx), peakIdx] = max(powerSpectrum);
    peakFrequency(stateIdx) = frequency(peakIdx);
end

[~, selectedIdx] = max(peakPower);
period = 1 / peakFrequency(selectedIdx);
end
