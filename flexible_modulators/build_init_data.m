clear
clc

%% Initial periodic orbit for the flexible modulator
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
dataDir = fullfile(scriptDir, 'data', 'common');
saveFile = fullfile(dataDir, 'initData.mat');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(dataDir);

Parameters = [1, 1];
y0 = [1; 0];

result = flexmod_forward_orbit(Parameters, y0, struct('systemName', 'base'));
if ~result.success
    error('flexmod_refactor:InitPeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        result.msg);
end

TS = {result.orbit.t, result.orbit.y};
period = result.features.period;
varAmp = reshape(result.features.state.amplitude, 1, []);
varMax = reshape(result.features.state.max, 1, []);
varMin = reshape(result.features.state.min, 1, []);

save(saveFile, 'TS', 'Parameters', 'period', 'varAmp', 'varMax', 'varMin');
fprintf('Saved init data: %s\n', saveFile);
