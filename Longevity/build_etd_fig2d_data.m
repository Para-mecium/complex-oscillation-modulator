clear
clc

%% File paths
scriptDir = fileparts(mfilename('fullpath'));
sdeDataFile = fullfile(scriptDir, 'sde_distribution.mat');
initDataFile = fullfile(scriptDir, 'initData.mat');
targetDataFile = fullfile(scriptDir, 'alpha_target_data.mat');

%% Load data
initData = load(initDataFile);
targetData = load(targetDataFile);
addpath(fullfile(scriptDir, 'SDE_simulation'));

%% Simulation settings
distributionSeeds = 1001:1100;
warmupT = 1000;
measureT = 60;
dt = 0.01;
sigma = 0.1;
noiseClass = repmat({'o'}, 1, 8);
measureCount = round(measureT / dt) + 1;

distributionOptions = struct();
distributionOptions.warmupT = warmupT;
distributionOptions.measureT = measureT;
distributionOptions.T = warmupT + measureT;
distributionOptions.dt = dt;
distributionOptions.sigma = sigma;
distributionOptions.noiseClass = noiseClass;

%% Prepare init state
initParams = reshape(double(initData.Parameters), 1, []);
initY0 = reshape(double(initData.TS{2}(1, :)), [], 1);

%% Prepare modulated state
modulatedParams = reshape(double(targetData.Parameters), 1, []);
modulatedY0 = reshape(double(targetData.TS{2}(1, :)), [], 1);

%% Simulate init distribution
initMinS = zeros(numel(distributionSeeds), 1);
initMinH = zeros(numel(distributionSeeds), 1);
initTrajectoryT = zeros(measureCount, 1);
initTrajectoryS = zeros(measureCount, numel(distributionSeeds));
initTrajectoryH = zeros(measureCount, numel(distributionSeeds));

disp('Start init simulation...')
for i = 1:numel(distributionSeeds)
    disp(['Simulating init seed ' num2str(i) '...'])
    distributionOptions.seed = distributionSeeds(i);
    initSim = Longevity_SDE(initParams, initY0, distributionOptions);
    initMask = initSim.t >= distributionOptions.warmupT;
    initMeasureT = initSim.t(initMask) - distributionOptions.warmupT;
    initMeasure = initSim.X(initMask, :);
    initTrajectoryT = initMeasureT;
    initTrajectoryS(:, i) = initMeasure(:, 3);
    initTrajectoryH(:, i) = initMeasure(:, 4);
    initMinS(i) = min(initMeasure(:, 3));
    initMinH(i) = min(initMeasure(:, 4));
end

init = struct();
init.distribution = struct();
init.distribution.minS = initMinS;
init.distribution.minH = initMinH;
init.trajectory = struct();
init.trajectory.t = initTrajectoryT;
init.trajectory.S = initTrajectoryS;
init.trajectory.H = initTrajectoryH;

%% Simulate modulated distribution
modulatedMinS = zeros(numel(distributionSeeds), 1);
modulatedMinH = zeros(numel(distributionSeeds), 1);
modulatedTrajectoryT = zeros(measureCount, 1);
modulatedTrajectoryS = zeros(measureCount, numel(distributionSeeds));
modulatedTrajectoryH = zeros(measureCount, numel(distributionSeeds));

disp('Start modulated simulation...')
for i = 1:numel(distributionSeeds)
    disp(['Simulating modulated seed ' num2str(i) '...'])
    distributionOptions.seed = distributionSeeds(i);
    modulatedSim = Longevity_SDE(modulatedParams, modulatedY0, distributionOptions);
    modulatedMask = modulatedSim.t >= distributionOptions.warmupT;
    modulatedMeasureT = modulatedSim.t(modulatedMask) - distributionOptions.warmupT;
    modulatedMeasure = modulatedSim.X(modulatedMask, :);
    modulatedTrajectoryT = modulatedMeasureT;
    modulatedTrajectoryS(:, i) = modulatedMeasure(:, 3);
    modulatedTrajectoryH(:, i) = modulatedMeasure(:, 4);
    modulatedMinS(i) = min(modulatedMeasure(:, 3));
    modulatedMinH(i) = min(modulatedMeasure(:, 4));
end

modulated = struct();
modulated.distribution = struct();
modulated.distribution.minS = modulatedMinS;
modulated.distribution.minH = modulatedMinH;
modulated.trajectory = struct();
modulated.trajectory.t = modulatedTrajectoryT;
modulated.trajectory.S = modulatedTrajectoryS;
modulated.trajectory.H = modulatedTrajectoryH;

%% Save data
save(sdeDataFile, 'init', 'modulated', 'distributionSeeds', 'distributionOptions');
