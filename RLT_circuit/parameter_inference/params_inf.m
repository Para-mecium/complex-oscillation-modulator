clear
clc
%% Data-driven procedure
% load data
load("initData_circuit.mat")

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

items_controlled = [1 4 5 6];
errBound = 1e-6;
continuationOptions = struct('initialLambdaStep', 0.01, 'predictorMode', 'constant');
Modtask = FMAM_ODE(sys,obs,StateView,items_per,items_controlled, [] ,errBound, ...
    'derivatives', derivatives, 'continuationOptions', continuationOptions);
% Modtask.fit()
Modtask.isPsiUpdated = true;
% Modtask.needLog = true;
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
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-6), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

poResult = extract_periodic_orbit(odeFunc, y0, Parameters, opts);

if ~poResult.has_orbit
    error('params_inf:PeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

State = state(obs, Parameters, poResult.orbit_t, poResult.orbit_y, M, PV);
State.updatePeriod();
State.updateVar2();

TS = {State.t, State.TS_var};
period = State.period;
varAmp = State.varAmp;
varMax = State.varMax;
varMin = State.varMin;

save(fullfile('RLT_circuit', 'learnedData_ODE.mat'), ...
    'TS', 'Parameters', 'period', 'varAmp', 'varMax', 'varMin');


%% 
figure
hold on
grid on
box on

State.TSplot('TS_var',1,'t')
State.TSplot('TS_var',2,'t')
State.TSplot('TS_var',3,'t')
