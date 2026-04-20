clear
clc

N = 100;
baseSeed = 0;
maxAttempts = 100;
sampledIdx = [1 4 5 6];
lowerBound = [500 0.5 0.5 0.5];
upperBound = [1500 2 2 2];

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
outputDir = fullfile(scriptDir, 'init_data_files');

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

oldFiles = dir(fullfile(outputDir, 'initData_*.mat'));
for i = 1:numel(oldFiles)
    delete(fullfile(oldFiles(i).folder, oldFiles(i).name));
end

data = load(fullfile(parameterInferenceDir, 'initData_ODE.mat'));
baseParameters = reshape(data.Parameters, 1, []);
y0 = data.TS{2}(1, :).';

acceptedCount = 0;

for attemptCount = 1:maxAttempts
    sampleSeed = baseSeed + attemptCount;
    rng(sampleSeed, 'twister');

    Parameters = baseParameters;
    Parameters(sampledIdx) = lowerBound + (upperBound - lowerBound) .* rand(1, numel(sampledIdx));

    result = circuit_forward_orbit(Parameters, y0, struct());
    if ~result.success
        continue
    end

    acceptedCount = acceptedCount + 1;

    TS = {result.orbit.t, result.orbit.y};
    seed = sampleSeed;
    period = result.features.period;
    varAmp = reshape(result.features.state.amplitude, 1, []);
    varMax = reshape(result.features.state.max, 1, []);
    varMin = reshape(result.features.state.min, 1, []);

    saveFile = fullfile(outputDir, sprintf('initData_%03d.mat', acceptedCount));
    save(saveFile, 'TS', 'Parameters', 'seed', 'period', 'varAmp', 'varMax', 'varMin');
    fprintf('Accepted %d/%d after %d attempts (seed=%d): %s\n', ...
        acceptedCount, N, attemptCount, sampleSeed, saveFile);

    if acceptedCount >= N
        break
    end
end

if acceptedCount < N
    error('build_init_data_files:InsufficientPeriodicOrbits', ...
        'Only found %d valid periodic orbits within %d attempts.', ...
        acceptedCount, maxAttempts);
end

fprintf('Saved %d valid periodic orbits after %d attempts.\n', ...
    acceptedCount, attemptCount);
