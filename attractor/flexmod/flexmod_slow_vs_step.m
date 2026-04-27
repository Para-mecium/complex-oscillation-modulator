% Base flexible modulator: compare slow parameter modulation and a sudden jump.

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(fileparts(scriptDir));
addpath(scriptDir, '-begin');
addpath(fullfile(repoDir, 'PO_extract'), '-begin');
clear forward_orbit build_model

%% 1. Local configuration
pathFile = fullfile(repoDir, 'flexible_modulators', 'data', 'fig3d', ...
    'curves', 'fig3d_iso_amplitude_A2p1.mat');
pathBranchName = 'rightBranch';

nConvergePeriods = 20;
nHoldPeriods = 5;
nRampPeriods = 50;
nAfterPeriods = 20;

RelTol = 1e-6;
AbsTol = 1e-9;
odeOptions = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
orbitOptions = struct('systemName', 'base');

obs = @(Y) mean(Y, 2);
obsLabel = 'Mean State';
stateLabels = {'X', 'Y'};
modelName = 'flexmod base';

%% 2. Load the modulation path endpoints
pathData = load(pathFile);
[p0, p1, pathLabel] = select_flexmod_path(pathData, pathBranchName);

%% 3. Extract the initial periodic orbit and converge on p0
fprintf('Extracting initial periodic orbit for %s...\n', modelName);
initialOrbit = forward_orbit(p0, [], orbitOptions);
if ~initialOrbit.success
    error('flexmod_slow_vs_step:InitialOrbitFailed', ...
        'Initial periodic-orbit extraction failed: %s', string(initialOrbit.message));
end

period0 = initialOrbit.orbit.period;
yCycle = initialOrbit.orbit.y(1, :).';
tConverge = nConvergePeriods * period0;

fprintf('Converging for %.3g time units before comparison...\n', tConverge);
[~, yConverge] = ode45(@(t, y) build_model(t, y, p0), ...
    [0, tConverge], yCycle, odeOptions);
yStart = yConverge(end, :).';

%% 4. Slow modulation and sudden jump
tHold = nHoldPeriods * period0;
tRampEnd = tHold + nRampPeriods * period0;
tTotal = tRampEnd + nAfterPeriods * period0;
tJump = tHold;

fprintf('Computing slow parameter modulation...\n');
[tSlow, ySlow] = ode45(@(t, y) build_model(t, y, ...
    linear_parameter_profile(t, p0, p1, tHold, tRampEnd)), ...
    [0, tTotal], yStart, odeOptions);

fprintf('Computing sudden parameter jump...\n');
[tStep, yStep] = ode45(@(t, y) build_model(t, y, ...
    step_parameter_profile(t, p0, p1, tJump)), ...
    [0, tTotal], yStart, odeOptions);

lambdaSlow = linear_lambda_profile(tSlow, tHold, tRampEnd);
lambdaStep = step_lambda_profile(tStep, tJump);
obsSlow = obs(ySlow);
obsStep = obs(yStep);

%% 5. Visualization
figure('Position', [50, 100, 1200, 600], ...
    'Name', sprintf('%s: Slow vs Sudden Parameter Change', modelName));

subplot(2, 2, 1);
plot(tSlow, obsSlow, 'b-', 'LineWidth', 1.2);
title('Slow Parameter Change', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(obsLabel, 'FontSize', 11);
grid on;
xlim([0, tTotal]);

subplot(2, 2, 3);
plot(tSlow, lambdaSlow, 'r-', 'LineWidth', 2);
xlabel('Time', 'FontSize', 11);
ylabel('\lambda(t)', 'FontSize', 11);
ylim([-0.05, 1.05]);
grid on;
xlim([0, tTotal]);

subplot(2, 2, 2);
plot(tStep, obsStep, 'b-', 'LineWidth', 1.2);
title('Sudden Parameter Jump', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(obsLabel, 'FontSize', 11);
grid on;
xlim([0, tTotal]);
xline(tJump, 'k--', 'Jump to p_1', ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 10);

subplot(2, 2, 4);
plot(tStep, lambdaStep, 'r-', 'LineWidth', 2);
xlabel('Time', 'FontSize', 11);
ylabel('\lambda(t)', 'FontSize', 11);
ylim([-0.05, 1.05]);
grid on;
xlim([0, tTotal]);

sgtitle(sprintf('%s (%s)', modelName, pathLabel), ...
    'FontSize', 16, 'FontWeight', 'bold');

plot_state_traces(tSlow, ySlow, tStep, yStep, stateLabels, modelName, tJump, tTotal);

%% ================= Local helpers ================= %%

function [p0, p1, pathLabel] = select_flexmod_path(pathData, branchName)
if ~isfield(pathData, 'seed') || ~isfield(pathData.seed, 'params') || ~isfield(pathData, branchName)
    error('select_flexmod_path:InvalidPathData', ...
        'Expected seed.params and branch %s in the selected path file.', branchName);
end

branch = pathData.(branchName);
if isempty(branch)
    error('select_flexmod_path:EmptyBranch', 'Selected branch is empty.');
end

p0 = reshape(pathData.seed.params, 1, []);
p1 = reshape(branch(end).params, 1, []);
pathLabel = sprintf('fig3d A2p1: seed -> %s(end)', branchName);
end

function p = linear_parameter_profile(t, p0, p1, tHold, tRampEnd)
lambda = linear_lambda_profile(t, tHold, tRampEnd);
p = (1 - lambda) * p0 + lambda * p1;
end

function p = step_parameter_profile(t, p0, p1, tJump)
if t < tJump
    p = p0;
else
    p = p1;
end
end

function lambda = linear_lambda_profile(t, tHold, tRampEnd)
lambda = zeros(size(t));
lambda(t >= tRampEnd) = 1;
mask = t >= tHold & t < tRampEnd;
lambda(mask) = (t(mask) - tHold) ./ (tRampEnd - tHold);
end

function lambda = step_lambda_profile(t, tJump)
lambda = double(t >= tJump);
end

function plot_state_traces(tSlow, ySlow, tStep, yStep, labels, modelName, tJump, tTotal)
figure('Position', [80, 120, 1200, 500], ...
    'Name', sprintf('%s: State Traces', modelName));
colors = lines(size(ySlow, 2));

subplot(1, 2, 1);
hold on;
for i = 1:size(ySlow, 2)
    plot(tSlow, ySlow(:, i), 'LineWidth', 1.2, 'Color', colors(i, :));
end
hold off;
title('Slow Parameter Change: State Variables', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time', 'FontSize', 11);
ylabel('State', 'FontSize', 11);
legend(labels, 'Location', 'best');
grid on;
xlim([0, tTotal]);

subplot(1, 2, 2);
hold on;
for i = 1:size(yStep, 2)
    plot(tStep, yStep(:, i), 'LineWidth', 1.2, 'Color', colors(i, :));
end
hold off;
title('Sudden Parameter Jump: State Variables', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time', 'FontSize', 11);
ylabel('State', 'FontSize', 11);
legend(labels, 'Location', 'best');
grid on;
xlim([0, tTotal]);
xline(tJump, 'k--', 'LineWidth', 1);
end
