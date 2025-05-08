clear
clc
%%
System
obs = [];

load("Longevity_data.mat")
ndim = 4;
ntst = 40;
iter = 9;

[t,T,TS_var] = xyzt(x,f,ndim,ntst,iter);
TS_var = TS_var';

params = zeros(1,13);
params(1) = 30.5; % alphaS
params(2) = 183; % alphaH
params(3) = 0.1; % alphaS0
params(4) = 0.1; % alphaH0
params(5) = x(end,iter); % beta1
params(6) = 3.7; % beta2
params(7) = 0.3; % deltam
params(8) = 0.2; % deltaS
params(9) = 3.8; % deltaH
params(10) = 326; % KH
params(11) = 185; % KS
params(12) = 3; % n1
params(13) = 4.8; % n2

M = 75; % truncation order
cd ../

% Set Primary variable
PV.name = 'var';
PV.idx = 1;
State = state(obs,params,t,TS_var,M,PV);

% Set modualtion target
varMin1_target = 300;
varMin2_target = 250;

items_per = struct;

items_per(1).prop = 'varMin';
items_per(1).idx = 3;
items_per(1).target = varMin1_target;

items_per(2).prop = 'varMin';
items_per(2).idx = 4;
items_per(2).target = varMin2_target;

items_controlled = [1 2];

% Set stepsize
stepsNum = 300;
max_stepsize_const = [1e-2 1e-2];
max_stepsize = max(abs([varMin1_target - State.varMin(1), varMin2_target - State.varMin(2)])/stepsNum,...
    max_stepsize_const);

errBound = 1e-6;
M = 75;
L = 500;

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

solution = Modtask.logs;
