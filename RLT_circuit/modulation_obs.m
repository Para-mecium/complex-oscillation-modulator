%% FMAM
% This block of codes performs a bi-modulation task of observable y = V1 + V2.
% We modulate the period, the amplitude of y through varying beta1 and beta2
clear

% load data
System
Observable

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

period_init = State.period;
obsAmp_init = State.obsAmp(1);

%%
% Set modualtion target
period_target = period_init;
obsAmp_target = 1.5;

items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period_target/(2*pi);

items_per(2).prop = 'obsAmp';
items_per(2).idx = 1;
items_per(2).target = obsAmp_target;

% Set control parameters
items_controlled = [1 4];
% Set stepsize
stepsNum = 200;
max_stepsize_const = [1e-3 1e-3];
max_stepsize = max(abs([(period_target-State.period)/(2*pi) ...
    obsAmp_target-State.obsAmp(1)])/stepsNum,...
    max_stepsize_const);

errBound = 1e-6;

Modtask = FMAM_ODE(sys,obs,State,items_per,items_controlled,max_stepsize,errBound);
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
State.updateObs2();

solution = Modtask.logs;