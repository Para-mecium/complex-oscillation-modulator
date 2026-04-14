plotsize=1;

figure()
for count=1:plotsize:size(solution,2)
%p_Psi=solution{4,count};
%q_Psi=solution{5,count};
p_variable=solution{6,count};
q_variable=solution{7,count};
phi = (0:L-1)'*2*pi/L;
[vc,vs] = Vec_CS(phi,M,L);
%Psi = vc*p_Psi+vs*q_Psi;
TS_variable = vc*p_variable+vs*q_variable;
plot(TS_variable(:,1),TS_variable(:,2))
%plot(phi,Psi)
hold on
end


%%
figure()
count=1;
%p_Psi=solution{4,count};
%q_Psi=solution{5,count};
p_variable=solution{6,count};
q_variable=solution{7,count};
phi = (0:L-1)'*2*pi/L;
[vc,vs] = Vec_CS(phi,M,L);
%Psi = vc*p_Psi+vs*q_Psi;
TS_variable = vc*p_variable+vs*q_variable;
plot(TS_variable(:,1),TS_variable(:,2))
hold on
count=size(solution,2);
%p_Psi=solution{4,count};
%q_Psi=solution{5,count};
p_variable=solution{6,count};
q_variable=solution{7,count};
phi = (0:L-1)'*2*pi/L;
[vc,vs] = Vec_CS(phi,M,L);
%Psi = vc*p_Psi+vs*q_Psi;
TS_variable = vc*p_variable+vs*q_variable;
plot(TS_variable(:,1),TS_variable(:,2))

%%
figure()
for count=1:plotsize:size(solution,2)
plot(solution{1,count}(Coe_Controlled{2,1}),solution{1,count}(Coe_Controlled{2,2}),'o')
hold on
end


%%
parameter=sol_para_AM(1,:);
[T,X] = solve_FHN(parameter);
figure()
plot(T,X)

function [T,X] = solve_FHN(parameter)
initial=[0,0];
tspan = 0:1:600;
abstol = 1e-2;
reltol = 1e-2;
ODEoptions = odeset('AbsTol',abstol,'RelTol',reltol);
[T,X] = ode45(@FHNode,tspan,initial,ODEoptions);

function output = FHNode(t,x)

    output = zeros(2,1);
    V=x(1);
    W=x(2);
    theta = parameter(1);
    gamma = parameter(2);
    I = parameter(3);
    epsilon = parameter(4);
    f11 = parameter(5);
    f12 = parameter(6);
    f21 = parameter(7);
    f22 = -parameter(5);
 
    output(1) = V.*(V-theta).*(1-V)-W+I+f11*(V-0.3068)+f12*(W-0.1227);
    output(2) = epsilon*(V-gamma*W)+f21*(V-0.3068)+f22*(W-0.1227);
end

end                                    