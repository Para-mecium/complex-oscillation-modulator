clear
clc
%% Data-driven procedure
%
needPath = true;
scriptDir = fileparts(mfilename('fullpath'));
% load data
load("initData_circuit.mat")

% Set modualtion target
items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period/(2*pi);

items_per(2).prop = 'varAmp';
items_per(2).idx = 1;
items_per(2).target = varAmp(1);

items_per(3).prop = 'varAmp';
items_per(3).idx = 2;
items_per(3).target = varAmp(2);

items_per(4).prop = 'varAmp';
items_per(4).idx = 3;
items_per(4).target = varAmp(3);

% Initialization
load("initData_ODE.mat")

cd ../
System 
obs = [];

params = Parameters;

t = TS{1};
TS_var = TS{2};
TS_obs = [];
M = 75; % truncation order
cd ../
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

% Set Primary variable
PV.name = 'var';
PV.idx = 1;
State = state(obs,params,t,TS_var,M,PV);
StateView = fmam_state_ops.solverViewFromState(State);

items_controlled = [1 4 5 6];
errBound = 1e-6;
continuationOptions = struct('initialLambdaStep', 0.01, 'predictorMode', 'constant');
Modtask = FMAM_ODE(sys,obs,StateView,items_per,items_controlled, [] ,errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
Modtask.psiUpdateMode = true;
Modtask.refreshPsiModeReferences();
Modtask.needLog = needPath;
%%
tic
Modtask.fit()
Modtask.step()
Modtask.errBound = 1e-12;
Modtask.fit()
Modtask.step()
elapsedTime = toc;  % 
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

StateView = Modtask.exportSolverView();
StateDerived = Modtask.exportDerivedView();

if needPath
    plotDataFile = fullfile(scriptDir, "params_modulation_path.mat");
    params_start_path = reshape(double(State.params), 1, []);
    params_end_path = reshape(double(StateView.params), 1, []);

    solution_path = Modtask.logs;
    if isempty(solution_path)
        error('params_inf:NoContinuationLogs', ...
            ['No continuation logs for params_modulation_path. Increase lambdaStepCap ' ...
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
odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter)];
searchWindow = max(10 * StateDerived.period, StateDerived.t(end));
poExtractDir = fullfile(pwd, 'PO_extract');
addpath(poExtractDir);

opts = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, searchWindow], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, opts);

if ~poResult.has_orbit
    error('params_inf:PeriodicOrbitGenerationFailed', ...
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

save(fullfile('RLT_circuit', 'learnedData_ODE.mat'), ...
    'TS', 'Parameters', 'period', 'varAmp', 'varMax', 'varMin');

%% 
figure
hold on
grid on
box on

plot(poResult.orbit_t, poResult.orbit_y(:, 1))
plot(poResult.orbit_t, poResult.orbit_y(:, 2))
plot(poResult.orbit_t, poResult.orbit_y(:, 3))
% plot(StateDerived.t, StateDerived.TS_var(:,1))
% plot(StateDerived.t, StateDerived.TS_var(:,2))
% plot(StateDerived.t, StateDerived.TS_var(:,3))