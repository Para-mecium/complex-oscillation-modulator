clear
clc

%% Initial periodic orbit for the circadian model
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
dataDir = fullfile(scriptDir, 'data', 'common');
saveFile = fullfile(dataDir, 'initData.mat');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
if ~isfolder(dataDir)
    mkdir(dataDir);
end

%% Initial parameter point and initial condition
Parameters = [1.0e-4, 0.05];
y0 = [0.10; 0.08; 0.07];

result = circadian_forward_orbit(Parameters, y0, struct());
if ~result.success
    error('circadian_refactor:InitPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        result.message);
end

%% Save the seed periodic orbit used by later fig. 5b scripts
TS = {result.orbit.t, result.orbit.y};
period = result.features.period;
obsAmp = reshape(result.features.observable.amplitude, 1, []);
obsMax = reshape(result.features.observable.max, 1, []);
obsMin = reshape(result.features.observable.min, 1, []);

save(saveFile, 'TS', 'Parameters', 'period', 'obsAmp', 'obsMax', 'obsMin');
fprintf('Saved init data: %s\n', saveFile);
