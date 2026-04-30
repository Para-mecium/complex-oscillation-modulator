clear
clc

%% User settings
settings = struct();
settings.randomSeed = 1;
settings.useParallel = true;
settings.numWorkers = 8;
settings.numRepeat = 100;
settings.numSourceNodeSequences = 5;
% Options: 'fixedNetwork' or 'independentNetworkDraws'
settings.repeatMode = 'fixedNetwork';
% Options: 'fixedAttempts' or 'untilSuccesses'
settings.repeatAccountingMode = 'fixedAttempts';
settings.maxAttemptsMultiplier = 10;
% Options: 'prefix', 'uniformRandomWithoutReplacement', 'uniformRandomWithReplacement'
settings.perturbedNodeSelectionMode = 'uniformRandomWithoutReplacement';
% Options: 'prefix', 'uniformRandomWithoutReplacement', 'uniformRandomWithReplacement'
settings.inputNodeSelectionMode = 'uniformRandomWithoutReplacement';
settings.weightPer = 0.3;
settings.outputSubdir = 'Ergodic data';
settings.solverName = 'ode15s';
settings.searchWindow = 3000;
% Options: any subset/permutation of {'SW', 'BA', 'ER'}
settings.nets = {'ER','BA','SW'};
settings.network = struct( ...
    'N', 50, ...
    'degAvg', 5, ...
    'K', 1, ...
    'baSeedDegree', 5, ...
    'swNeighborCount', 2, ...
    'swRewireProb', 0.15);
settings.modelParams = struct( ...
    'theta', 0.5, ...
    'gamma', 1, ...
    'epsilon', 1e-2, ...
    'tau', 1, ...
    'I', 0.38);

modelSpec = create_model_spec(settings.modelParams, settings.searchWindow);

tic;
networkexp.run_ergodic_experiment(settings, modelSpec);
elapsedTime = toc;
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

function modelSpec = create_model_spec(modelParams, searchWindow)
modelSpec = struct( ...
    'name', 'FHN', ...
    'makeInitialState', @(N) 0.21 * ones(2 * N, 1), ...
    'buildParams', @(omega) { ...
        modelParams.theta, ...
        modelParams.gamma, ...
        modelParams.epsilon, ...
        modelParams.tau, ...
        modelParams.I, ...
        sparse(omega)}, ...
    'buildSolverOptions', @(params) networkexp.build_fhn_orbit_detection_options( ...
        networkexp.build_fhn_solver_options(params), searchWindow), ...
    'rhs', @fhn_rhs);
end

function dydt = fhn_rhs(~, x, params)
theta = params{1};
gamma = params{2};
epsilon = params{3};
tau = params{4};
I = params{5};
matAdj = sparse(params{6});

V = 0;
N = size(matAdj, 1);
dydt = zeros(2 * N, 1);
v = x(1:N);
w = x(N + 1:end);
sigma = 1 ./ (1 + exp(-tau * (v - V)));

dydt(1:N) = v .* (v - theta) .* (1 - v) - w + I + matAdj * sigma;
dydt(N + 1:end) = epsilon * (v - gamma * w);
end
