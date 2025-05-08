%% FMAM
% This block of codes performs a bi-modulation task.
% We modulate the period, the amplitude of V1 through varying beta1 and 
% beta2.
clear

% load data
System
obs = [];

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
varAmp1_init = State.varAmp(1);
%%
% Set modualtion target
period_target = period_init*2;
varAmp1_target = varAmp1_init;

items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period_target/(2*pi);

items_per(2).prop = 'p_var';
items_per(2).idx = [2,1];
items_per(2).target = varAmp1_target;

% items_per(3).prop = 'varAmp';
% items_per(3).idx = 2;
% items_per(3).target = State.params(3)*25;

% Set control parameters
items_controlled = [1 4];
% Set stepsize
stepsNum = 100;
max_stepsize_const = [1e-3 1e-2];
max_stepsize = max(abs([(period_target-State.period)/(2*pi) ...
    varAmp1_target-State.varAmp(1)])/stepsNum,...
    max_stepsize_const);

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

solution = Modtask.logs;