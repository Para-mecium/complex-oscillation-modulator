clear
clc
needPath = true;

%% File paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'initData.mat');
saveFile = fullfile(scriptDir, 'beta_target_data.mat');
plotDataFile = fullfile(scriptDir, 'beta_modulation_path.mat');

%% Load initialization data
obs = [];
M = 75;
PV.name = 'var';
PV.idx = 1;

addpath(fullfile(repoDir, 'PO_extract'));

initData = load(initDataFile);
params = initData.Parameters;
t = initData.TS{1};
TS_var = initData.TS{2};

%% Build FMAM inputs
targetMinSH = [300 300];
controlledIdx = [5 6];
continuationOptions = struct('initialLambdaStep', 0.01, 'predictorMode', 'auto');

sys = build_longevity_system();
derivatives = build_symbolic_derivatives(sys, obs, numel(params));
State = state(obs, params, t, TS_var, M, PV);
StateView = fmam_state_ops.solverViewFromState(State);

items_per = struct;
items_per(1).prop = 'varMin';
items_per(1).idx = 3;
items_per(1).target = targetMinSH(1);

items_per(2).prop = 'varMin';
items_per(2).idx = 4;
items_per(2).target = targetMinSH(2);

items_controlled = controlledIdx;
errBound = 1e-6;
Modtask = FMAM_ODE(sys, obs, StateView, items_per, items_controlled, [], errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
Modtask.psiUpdateMode = true;
Modtask.refreshPsiModeReferences();
Modtask.needLog = needPath;

%% Run FMAM continuation
tic
Modtask.fit()
Modtask.step()
Modtask.errBound = 1e-12;
Modtask.fit()

elapsedTime = toc;
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

StateView = Modtask.exportSolverView();
StateDerived = Modtask.exportDerivedView();

%% Save continuation path
if needPath
    params_start_path = reshape(double(State.params), 1, []);
    params_end_path = reshape(double(StateView.params), 1, []);

    solution_path = Modtask.logs;
    if isempty(solution_path)
        error('build_figS13a_cache_fmam:NoContinuationLogs', ...
            ['No continuation logs for figS13a_params_modulation_path. Increase lambdaStepCap ' ...
             'or verify the continuation target differs from the initial state.']);
    end

    curve_params_path = zeros(numel(solution_path), numel(solution_path(1).params));
    for i = 1:numel(solution_path)
        curve_params_path(i, :) = reshape(double(solution_path(i).params), 1, []);
    end

    pathData = struct();
    pathData.curve_params = curve_params_path;
    pathData.params_start = params_start_path;
    pathData.params_end = params_end_path;
    save(plotDataFile, '-struct', 'pathData');
    fprintf('Saved continuation path cache: %s\n', plotDataFile);
end

%% Extract ODE periodic orbit at inferred parameters
Parameters = StateView.params;
inferredTS = StateDerived.TS_var;
y0 = inferredTS(1, :).';
searchWindow = max(10 * StateDerived.period, StateDerived.t(end));
orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 400], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);
opts = orbitOptions;
opts.tspan = [0, max(searchWindow, orbitOptions.tspan(end))];
odeFunc = @(t, y, parameter) ode_rhs_from_sys(sys, y, parameter);

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, opts);
if ~poResult.has_orbit
    error('build_figS13a_cache_fmam:PeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

TS = {poResult.orbit_t, poResult.orbit_y};
period = poFeatures.period;
varAmp = reshape(poFeatures.state.amplitude, 1, []);
varMax = reshape(poFeatures.state.max, 1, []);
varMin = reshape(poFeatures.state.min, 1, []);

save(saveFile, 'TS', 'Parameters', 'period', 'varAmp', 'varMax', 'varMin');
fprintf('Saved learned data: %s\n', saveFile);

%%
function sys = build_longevity_system()
sys = cell(1, 4);
sys{1} = @dmSdt;
sys{2} = @dmHdt;
sys{3} = @dSdt;
sys{4} = @dHdt;
end

function dydt = ode_rhs_from_sys(sys, y, parameter)
yRow = reshape(y, 1, []);
dydt = zeros(numel(sys), 1);
for i = 1:numel(sys)
    dydt(i) = sys{i}(yRow, parameter);
end
end

function output = dmSdt(TS_variable, parameter)
alphaS = parameter(1);
alphaS0 = parameter(3);
deltam = parameter(7);
KH = parameter(10);
n1 = parameter(12);

mS = TS_variable(:, 1);
H = TS_variable(:, 4);

output = alphaS0 + alphaS * H.^n1 ./ (KH^n1 + H.^n1) - deltam * mS;
end

function output = dmHdt(TS_variable, parameter)
alphaH = parameter(2);
alphaH0 = parameter(4);
deltam = parameter(7);
KS = parameter(11);
n2 = parameter(13);

mH = TS_variable(:, 2);
S = TS_variable(:, 3);

output = alphaH0 + alphaH * KS^n2 ./ (KS^n2 + S.^n2) - deltam * mH;
end

function output = dSdt(TS_variable, parameter)
betaS = parameter(5);
deltaS = parameter(8);

mS = TS_variable(:, 1);
S = TS_variable(:, 3);

output = betaS * mS - deltaS * S;
end

function output = dHdt(TS_variable, parameter)
betaH = parameter(6);
deltaH = parameter(9);

mH = TS_variable(:, 2);
H = TS_variable(:, 4);

output = betaH * mH - deltaH * H;
end
