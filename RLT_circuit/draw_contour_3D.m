lineWidth = 0.5;
mksize = 10;
stepsize = 10;
%% ODE plot 3-D
% plot 3-params curve
load("contour T = T0-6 ~ T0+25, Ax1 = A10, Ax2 = A20, Imax = 2.97mA.mat")
curve = solution;

curve_params = zeros(size(curve,2),size(curve(1).params,2));

for i = 1:size(curve,2)
    curve_params(i,:) = curve(i).params;
end

figure

plot3(1./curve_params(1:stepsize:end,4),1./curve_params(1:stepsize:end,5),curve_params(1:stepsize:end,1),LineWidth=lineWidth*2);
grid on
hold on
box on

% load("contour data -> T = 10, Ax1 = 10, Ax2 = 25.mat")
% curve = solution;
% 
% curve_params = zeros(size(curve,2),size(curve(1).params,2));
% for i = 1:size(curve,2)
%     curve_params(i,:) = curve(i).params;
% end
% plot3(1./curve_params(:,4),1./curve_params(:,5),curve_params(:,1),LineWidth=lineWidth*2);

params_plot = {};
params_plot{1} = [946.2608438	6.6	0.06	2.247628785	0.405358612	1.118850063	0.00297];
params_plot{2} = [1036.49162720506	6.60000000000000	0.0600000000000000	1/0.273726121820823	1/1.78970815080701	1/0.893774807719991	0.00297000000000000];
params_plot{3} = [1195.33831135089	6.60000000000000	0.0600000000000000	1/0.219022230994216	1/1.69303561778188	1/0.893774807719991	0.00297000000000000];
% plot params points for T = 10, 15, 20
for i = 1:3
    para = params_plot{i};
    plot3(ones(1,2)*para(4),ones(1,2)*para(5),[para(1) -1],MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
    plot3(ones(1,2)*para(4),[para(5) 2],ones(1,2)*para(1),MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
end

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

%% plot 2-params curve (Rc, C3)
load("contour T = T0-6 ~ T0+25, Ax1 = A10, Ax2 = A20, Imax = 2.97mA.mat")
curve = solution;

curve_params = zeros(size(curve,2),size(curve(1).params,2));

for i = 1:size(curve,2)
    curve_params(i,:) = curve(i).params;
end

figure

plot(1./curve_params(1:stepsize:end,6),curve_params(1:stepsize:end,1),LineWidth=lineWidth*2)

grid on
box on
hold on

params_plot = {};
params_plot{1} = [946.2608438	6.6	0.06	2.247628785	0.405358612	1.118850063	0.00297];
params_plot{2} = [1036.49162720506	6.60000000000000	0.0600000000000000	1/0.273726121820823	1/1.78970815080701	1/0.893774807719991	0.00297000000000000];
params_plot{3} = [1195.33831135089	6.60000000000000	0.0600000000000000	1/0.219022230994216	1/1.69303561778188	1/0.893774807719991	0.00297000000000000];
% plot params points for T = 10, 15, 20
for i = 1:3
    para = params_plot{i};
    plot(para(6),para(1),MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
end
% plot(1.118850063,946.2608438,MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)

axis([1 1.2 900 1500])

xlabel('{\itC}_3 (\muF)','Fontname','Arial')
ylabel('{\itR_C} (k\Omega)','Fontname','Arial')

XTick_2D = 1:0.05:1.2;
YTick_2D = 0:300:1500;

XTickLabel_2D = {'1.00','1,05','1.10','1.15','1.20'};
YTickLabel_2D = {'0','300','600','900','1200','1500'};

set(gca,'FontSize',10,'XTick',XTick_2D,'YTick',YTick_2D,'XTickLabel',XTickLabel_2D,'YTickLabel',YTickLabel_2D)
% view(-65, 40); 