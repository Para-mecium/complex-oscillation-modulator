clear
clc

%% File paths
scriptDir = fileparts(mfilename('fullpath'));
sdeDataFile = fullfile(scriptDir, 'sde_representative.mat');
initDataFile = fullfile(scriptDir, 'initData.mat');
targetDataFile = fullfile(scriptDir, 'alpha_target_data.mat');

%% Load data
initData = load(initDataFile);
targetData = load(targetDataFile);
addpath(fullfile(scriptDir, 'SDE_simulation'));

%% Simulation settings
seed = 1;
warmupT = 1000;
measureT = 60;

Options = struct();
Options.seed = seed;
Options.warmupT = warmupT;
Options.measureT = measureT;
Options.T = warmupT + measureT;
Options.dt = 0.001;
Options.sigma = 0.1;
Options.noiseClass = repmat({'o'}, 1, 8);

%% Prepare init state
initParams = reshape(double(initData.Parameters), 1, []);
initY0 = reshape(double(initData.TS{2}(1, :)), [], 1);

%% Prepare modulated state
modulatedParams = reshape(double(targetData.Parameters), 1, []);
modulatedY0 = reshape(double(targetData.TS{2}(1, :)), [], 1);

%% Simulate init representative trace
initSim = Longevity_SDE(initParams, initY0, Options);
initMask = initSim.t >= Options.warmupT;

init = struct();
init.representative = struct();
init.representative.t = initSim.t(initMask) - Options.warmupT;
init.representative.X = initSim.X(initMask, :);

%% Simulate modulated representative trace
modulatedSim = Longevity_SDE(modulatedParams, modulatedY0, Options);
modulatedMask = modulatedSim.t >= Options.warmupT;

modulated = struct();
modulated.representative = struct();
modulated.representative.t = modulatedSim.t(modulatedMask) - Options.warmupT;
modulated.representative.X = modulatedSim.X(modulatedMask, :);

%% Save data
save(sdeDataFile, 'init', 'modulated', 'seed', 'Options');
