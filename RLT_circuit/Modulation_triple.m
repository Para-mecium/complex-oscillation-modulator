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
period_init = State.period;
varAmp1_init = State.varAmp(1);
varAmp2_init = State.varAmp(2);
%%
% Set modualtion target
period_target = period_init*2;
varAmp1_target = varAmp1_init;
varAmp2_target = varAmp2_init;
items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period_target/(2*pi);

items_per(2).prop = 'p_var';
items_per(2).idx = [2,1];
items_per(2).target = varAmp1_target;

items_per(3).prop = 'varAmp';
items_per(3).idx = 2;
items_per(3).target = varAmp2_target;

% Set control parameters
items_controlled = [1 4 5];
% Set stepsize
stepsNum = 100;
max_stepsize_const = [1e-3 1e-2 1e-2];
max_stepsize = max(abs([(period_target-State.period)/(2*pi) ...
    varAmp1_target-State.varAmp(1) varAmp2_target-State.varAmp(2)])/stepsNum,...
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