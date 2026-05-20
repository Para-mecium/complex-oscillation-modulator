clear
clc

N = 100;
maxAttempts = 100;
scaleLevels = [1 0.5 0.75 1.25 1.5];
baseSeeds = [0 1000 2000 3000 4000];
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

data = load(fullfile(parameterInferenceDir, 'initData_ODE.mat'));
baseParameters = reshape(data.Parameters, 1, []);
y0 = data.TS{2}(1, :).';

Width = upperBound - lowerBound;

for scaleIdx = 1:numel(scaleLevels)
    scale = scaleLevels(scaleIdx);
    baseSeed = baseSeeds(scaleIdx);
    scaledLowerBound = lowerBound;
    scaledUpperBound = lowerBound + scale * Width;

    scaleDir = fullfile(outputDir, scale_dir_name(scale));
    if ~isfolder(scaleDir)
        mkdir(scaleDir);
    end

    oldFiles = dir(fullfile(scaleDir, 'initData_*.mat'));
    for i = 1:numel(oldFiles)
        delete(fullfile(oldFiles(i).folder, oldFiles(i).name));
    end

    acceptedCount = 0;

    fprintf('Generating initial data for scale=%.6g, seedBase=%d, region=[%s] to [%s]\n', ...
        scale, baseSeed, num2str(scaledLowerBound), num2str(scaledUpperBound));

    for attemptCount = 1:maxAttempts
        sampleSeed = baseSeed + attemptCount;
        rng(sampleSeed, 'twister');

        Parameters = baseParameters;
        Parameters(sampledIdx) = scaledLowerBound + ...
            (scaledUpperBound - scaledLowerBound) .* rand(1, numel(sampledIdx));

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

        saveFile = fullfile(scaleDir, sprintf('initData_%03d.mat', acceptedCount));
        save(saveFile, 'TS', 'Parameters', 'seed', 'scale', ...
            'period', 'varAmp', 'varMax', 'varMin');
        fprintf('Accepted %d/%d after %d attempts (seed=%d): %s\n', ...
            acceptedCount, N, attemptCount, sampleSeed, saveFile);

        if acceptedCount >= N
            break
        end
    end

    if acceptedCount < N
        error('build_init_data_files:InsufficientPeriodicOrbits', ...
            'Scale %.6g only found %d valid periodic orbits within %d attempts.', ...
            scale, acceptedCount, maxAttempts);
    end

    fprintf('Saved %d valid periodic orbits for scale=%.6g after %d attempts.\n', ...
        acceptedCount, scale, attemptCount);
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end
