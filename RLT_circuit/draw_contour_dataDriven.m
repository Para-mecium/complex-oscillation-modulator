lineWidth = 0.5;
mksize = 10;
stepsize = 10;
colors = {'#1E4C9C', '#345D82', '#3371B3', '#AED4E5'};
color_l = colors{1};
%% plot 3-params curve (Rc, C1, C2)
load("contour dataDriven.mat")
curve = solution;

curve_params = zeros(size(curve,2),size(curve(1).params,2));

for i = 1:691
    curve_params(i,:) = curve(i).params;
end

curve_params = [curve_params; curve(end).params];

figure

plot3(1./curve_params(1:stepsize:end,4),1./curve_params(1:stepsize:end,5),curve_params(1:stepsize:end,1),LineWidth=lineWidth*2)
grid on
hold on
box on

% plot params points for T = 10, 15, 20

% T = 10 A_1 = 10 A_2 = 25
% 1001.25391586687	6.60000000000000	0.0600000000000000	0.527083357162762	2.97348864805823	0.951756736096122	0.00300000000000

% T = 15 A_1 = 10 A_2 = 25
% 1069.86029432421	6.60000000000000	0.0600000000000000	0.308371740150248	2.05270354797472	0.951756736096122	0.00300000000000000

% T = 20 A_1 = 10 A_2 = 25
% 1210.98355396663	6.60000000000000	0.0600000000000000	0.241942639153974	1.86970852879969	0.951756736096122	0.00300000000000000

plot3(ones(1,2)*2.247628785,ones(1,2)*0.405358612,[946.2608438 -1],MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
plot3([2.247628785 7],ones(1,2)*0.405358612,ones(1,2)*946.2608438,MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
plot3(ones(1,2)*2.247628785,[0.405358612 2],ones(1,2)*946.2608438,MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
plot3(1,1,500,MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)

axis([1 6 0 1 0 1500])

xlabel('{\itC}_1 (\muF)','Fontname','Arial')
ylabel('{\itC}_2 (\muF)','Fontname','Arial')
zlabel('{\itR_C} (k\Omega)','Fontname','Arial')

XTick_3D = 0:1:6;
YTick_3D = 0:0.2:1;
ZTick_3D = 0:300:1500;
XTickLabel_3D = {'0','1','2','3','4','5','6'};
YTickLabel_3D = {'0','0.2','0.4','0.6','0.8','1.0'};
ZTickLabel_3D = {'0','300','600','900','1200','1500'};

set(gca,'FontSize',10,'XTick',XTick_3D,'YTick',YTick_3D,'ZTick',ZTick_3D,'XTickLabel',XTickLabel_3D,'YTickLabel',YTickLabel_3D,'ZTickLabel',ZTickLabel_3D)
view(-70, 40); 

%% %% plot 2-params curve (Rc, C3)
load("contour dataDriven.mat")
curve = solution;

curve_params = zeros(size(curve,2),size(curve(1).params,2));

for i = 1:691
    curve_params(i,:) = curve(i).params;
end
curve_params = [[500 6.6 0.06 1 1 1 2.97e-3];curve_params; curve(end).params];

figure

plot(1./curve_params(1:stepsize:end,6),curve_params(1:stepsize:end,1),LineWidth=lineWidth*2)
grid on
box on
hold on

plot(1/curve_params(end,6),curve_params(end,1),MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
plot(1,500,MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)

axis([1 1.2 0 1500])

xlabel('{\itC}_3 (\muF)','Fontname','Arial')
ylabel('{\itR_C} (k\Omega)','Fontname','Arial')

XTick_2D = 1:0.05:1.2;
YTick_2D = 0:300:1500;

XTickLabel_2D = {'1.00','1,05','1.10','1.15','1.20'};
YTickLabel_2D = {'0','300','600','900','1200','1500'};

set(gca,'FontSize',10,'XTick',XTick_2D,'YTick',YTick_2D,'XTickLabel',XTickLabel_2D,'YTickLabel',YTickLabel_2D)
% view(-65, 40); 