
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


PV = cell(2,1);
PV{1,1} = 'state';
PV{2,1} = 1;

%%
value_target = [0.088907389856442;0.028086062518803;80];
maxstepsize = [0.1;0.1;0.1];
ExpectedError = 1e-6;
M = 200;
L = 500;

[solution,~]=FMAM_ODE(sys,obs,TimeSeries,Parameters,...
    Coe_Controlled,Coe_Perturb,value_target,maxstepsize,ExpectedError,M,L,PV);

count1=65;
for i=1:count1
Parameters_FM=cell2mat(solution(1,end));
sol_para_FM(i,:)=cell2mat(solution(1,end));
FourierCoeffs=solution(2:7,end);
FourierCoeffs{7}=[];
FourierCoeffs{8}=[];
value_target = [0.088907389856442;0.028086062518803;80-i];
[solution,~]=FMAM_ODE_FC(sys,obs,FourierCoeffs,Parameters_FM,Coe_Controlled,Coe_Perturb,value_target,maxstepsize,ExpectedError,M,L,PV);
end

%%
value_target = [0.088907389856442*0.7;0.028086062518803*0.7;32.6661806687750];
maxstepsize = [0.001;0.001;0.001];
ExpectedError = 1e-6;
M = 200;
L = 500;
[solution,~]=FMAM_ODE(sys,obs,TimeSeries,Parameters,...
    Coe_Controlled,Coe_Perturb,value_target,maxstepsize,ExpectedError,M,L,PV);

count2=80;
for j=1:count2
Parameters_AM=cell2mat(solution(1,end));
sol_para_AM(j,:)=cell2mat(solution(1,end));
FourierCoeffs=solution(2:7,end);
FourierCoeffs{7}=[];
FourierCoeffs{8}=[];
k=0.7+0.01*j;
value_target = [0.088907389856442*k;0.028086062518803*k;32.6661806687750];
[solution,~]=FMAM_ODE_FC(sys,obs,FourierCoeffs,Parameters_AM,Coe_Controlled,Coe_Perturb,value_target,maxstepsize,ExpectedError,M,L,PV);
end
