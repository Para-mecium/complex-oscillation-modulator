
System
obs = [];

load("FHN.mat")
ndim = 2;
ntst = 40;
iter = 12;

cd ../

[t,T,xx] = xyzt(x,f,ndim,ntst,iter);
TimeSeries = [t,xx'];

Parameters = zeros(1,7);
Parameters(1) = 0.2; % theta
Parameters(2) = 2.5; % gamma
Parameters(3) = 0.1; % I
Parameters(4) = 0.1; % epsilon
Parameters(5) = 0; % f11
Parameters(6) = 0; % f12
Parameters(7) = 0; % f21


Coe_Controlled = cell(2,3);
Coe_Controlled{1,1} = 'parameter';
Coe_Controlled{1,2} = 'parameter';
Coe_Controlled{1,3} = 'parameter';
Coe_Controlled{2,1} = 5; % f11
Coe_Controlled{2,2} = 6; % f12
Coe_Controlled{2,3} = 7; % f21

Coe_Perturb = cell(3,3);
Coe_Perturb{1,1} = 'state';
Coe_Perturb{2,1} = 1; % V
Coe_Perturb{3,1} = 1; % amplitude

Coe_Perturb{1,2} = 'state';
Coe_Perturb{2,2} = 2; % W
Coe_Perturb{3,2} = 1; % amplitude

Coe_Perturb{1,3} = 'parameter';
Coe_Perturb{2,3} = 0; % period
Coe_Perturb{3,3} = []; 

value_target = [0.088907389856442;0.028086062518803;32.6661806687750/2];
maxstepsize = [0.1;0.1;0.1];
ExpectedError = 1e-6;
M = 240;
L = 500;

PV = cell(2,1);
PV{1,1} = 'state';
PV{2,1} = 1;

[solution,Error]=FMAM_ODE(sys,obs,TimeSeries,Parameters,...
    Coe_Controlled,Coe_Perturb,value_target,maxstepsize,ExpectedError,M,L,PV);
