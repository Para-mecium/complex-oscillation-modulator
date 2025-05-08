clear
clc
%% Data-driven procedure
% load data
load("initData_circuit.mat")

% Initialization
System 
obs = [];
load("initData_ODE.mat")

params = Parameters;

t = TS{1};
TS_var = TS{2};
TS_obs = [];
M = 75; % truncation order
cd ../

% Set Primary variable
PV.name = 'var';
PV.idx = 1;
State = state(obs,params,t,TS_var,M,PV);

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
max_stepsize = abs([(period - State.period)/(2*pi) State.varAmp - varAmp])/700;
errBound = 1e-6;
Modtask = FMAM_ODE(sys,obs,State,items_per,items_controlled,max_stepsize,errBound);
% Modtask.fit()
Modtask.isPsiUpdated = true;
Modtask.fit()
Modtask.needLog = true;
%%
tic
Modtask.fit()
Modtask.step()
elapsedTime = toc;  % 
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

State.updatePeriod();
State.updateVar2();

solution = Modtask.logs;

%
% items_per(4).prop = 'varAmp';
% items_per(4).idx = 3;
% items_per(4).target = varAmp(3)*params(3);
% 
% items_controlled = [1 4 5 6];
% max_stepsize = [1e-2 1e-3 1e-3 1e-3];
% errBound = 1e-6;
% Modtask = FMAM_ODE(sys,[],State,items_per,items_controlled,max_stepsize,errBound);
% 
% %
% tic
% Modtask.fit()
% Modtask.step()
% elapsedTime = toc;  % 
% disp(['Computing time: ', num2str(elapsedTime), ' seconds']);
% 
% State.updatePeriod();
% State.updateVar2();

