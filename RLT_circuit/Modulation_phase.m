%% FMAM
% This block of codes performs a tri-modulation task.
% We modulate the period, the amplitude of V1 and the amplitude of V2
% through varying Rc, beta1 and beta2
clear

% load data
System
obs = [];

% load("T_A1_A2_A3.mat")
% load("initData (ODE, Imax = 2.97mA).mat")
load("learnedData_ODE.mat")

params = Parameters;

t = TS{1};
TS_var = TS{2};
TS_obs = [];
M = 75;
cd ../

% Set Primary variable
PV.name = 'var';
PV.idx = 1;
State = state(obs,params,t,TS_var,M,PV);
%%
% Set modualtion target
period_target = 30;
amplitude1_target = State.params(3)*10;
items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period_target/(2*pi);

items_per(2).prop = 'p_var';
items_per(2).idx = [2,1];
items_per(2).target = amplitude1_target;

% Set control parameters
items_controlled = [1 4];
% Set stepsize
max_stepsize = abs([(period_target-State.period)/(2*pi) amplitude1_target-State.varAmp(1)])/100;
errBound = 1e-6;

Modtask = FMAM_ODE(sys,[],State,items_per,items_controlled,max_stepsize,errBound);
Modtask.isPsiUpdated = true;
Modtask.needLog = true;
%%
tic
Modtask.fit()
Modtask.step()
elapsedTime = toc;  % 
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

State.updatePeriod();
State.updateVar2();

%%
% Set modualtion target
phase12_target = 25;
items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period_target/(2*pi);

items_per(2).prop = 'varPhase';
items_per(2).idx = [1 2];
items_per(2).target = phase12_target;

% Set control parameters
items_controlled = [4 5];
% Set stepsize
% max_stepsize = abs([(period_target-State.period)/(2*pi) phase12_target-State.varPhase(1,2)])/200;
max_stepsize = [0.1 0.1];
errBound = 1e-6;

Modtask = FMAM_ODE(sys,[],State,items_per,items_controlled,max_stepsize,errBound);
Modtask.isPsiUpdated = true;
Modtask.needLog = true;

tic
Modtask.fit()
Modtask.step()
elapsedTime = toc;  % 
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

State.updatePeriod();
State.updateVar2();