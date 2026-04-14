System
obs = [];

load("cancer.mat")
ndim = 2;
ntst = 40;
iter = 10;

cd ../

[t,T,xx] = xyzt(x,f,ndim,ntst,iter);
TimeSeries = [t,xx'];

Parameters = zeros(1,8);
Parameters(1) = 0.35; % alpha
Parameters(2) = 5; % kappa
Parameters(3) = 1; % gamma1
Parameters(4) = 2.5; % gamma2
Parameters(5) = 0.08;% epsilon
Parameters(6) = 0; % f11
Parameters(7) = 0; % f12
Parameters(8) = 0; % f21


Coe_Controlled = cell(2,3);
Coe_Controlled{1,1} = 'parameter';
Coe_Controlled{1,2} = 'parameter';
Coe_Controlled{1,3} = 'parameter';
Coe_Controlled{2,1} = 6; % f11
Coe_Controlled{2,2} = 7; % f12
Coe_Controlled{2,3} = 8; % f21

Coe_Perturb = cell(3,3);
Coe_Perturb{1,1} = 'state';
Coe_Perturb{2,1} = 1; % x
Coe_Perturb{3,1} = 1; % amplitude

Coe_Perturb{1,2} = 'state';
Coe_Perturb{2,2} = 2; % y
Coe_Perturb{3,2} = 1; % amplitude

Coe_Perturb{1,3} = 'parameter';
Coe_Perturb{2,3} = 0; % period
Coe_Perturb{3,3} = []; 

value_target = [1.218533129973287;0.652707264091526;10];%4.308975306457719
maxstepsize = [0.001;0.001;0.001];
ExpectedError = 1e-6;
M = 200;
L = 500;

PV = cell(2,1);
PV{1,1} = 'state';
PV{2,1} = 1;

[solution,Error]=FMAM_ODE(sys,obs,TimeSeries,Parameters,...
    Coe_Controlled,Coe_Perturb,value_target,maxstepsize,ExpectedError,M,L,PV);