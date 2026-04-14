cd ../
plotsize=1;

%% Longevity Phase diagram 1
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
plot(TS_variable(:,3),TS_variable(:,4))
%plot(phi,Psi)
hold on
end
box on
grid on
xlabel('S','Fontname','Arial')
ylabel('H','Fontname','Arial')
title('Longevity Phase diagram')
set(gca,'Fontsize',18)

%% Longevity Phase diagram 2
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
plot([TS_variable(:,3);TS_variable(1,3)],[TS_variable(:,4);TS_variable(1,4)],'LineWidth',2)
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
plot([TS_variable(:,3);TS_variable(1,3)],[TS_variable(:,4);TS_variable(1,4)],'LineWidth',2)
box on
grid on
xlabel('S','Fontname','Arial')
ylabel('H','Fontname','Arial')
title('Longevity Phase Diagram')
legend('before','after')
set(gca,'Fontsize',18)

%% Parameters Diagram (alphaS versus alphaH)
figure()
for count=1:plotsize:size(solution,2)
plot(solution{1,count}(Coe_Controlled{2,1}),solution{1,count}(Coe_Controlled{2,2}),'o')
hold on
end
grid on
set(gca,'Fontsize',18)
title('Parameters Diagram')
xlabel('$\alpha_S$','Fontname','Arial','Interpreter','latex')
ylabel('$\alpha_H$','Fontname','Arial','Interpreter','latex')

