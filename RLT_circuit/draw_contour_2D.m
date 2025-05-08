lineWidth = 0.5;
mksize = 10;
stepsize = 1;

%% plot 3-params curve (Rc, C1, C2)
load("contour T = T0-7 ~ T0+25, A1 = A10.mat")
curve = solution;

curve_params = zeros(size(curve,2),size(curve(1).params,2));

for i = 1:size(curve,2)
    curve_params(i,:) = curve(i).params;
end

figure
plot3(1./curve_params(1:stepsize:end,4),1/curve_params(1,5)*ones(size(curve(1:stepsize:end),2),1),curve_params(1:stepsize:end,1),LineWidth=2*lineWidth)

grid on
hold on 
box on

surf([-1 7;-1 7],0.405358612*ones(2,2),[0 0;1500 1500],facecolor=[0.5 0.5 0.5],facealpha=0.25,EdgeColor='none')

params_plot = {};
params_plot{1} = [946.2608438	6.6	0.06	2.247628785	0.405358612	1.118850063	0.00297];
params_plot{2} = [1180.85370223743	6.60000000000000	0.0600000000000000	1/0.323329693749483	1/2.46695141192899	1/0.893774807719991	0.00297000000000000];
params_plot{3} = [1389.62876468867	6.60000000000000	0.0600000000000000	1/0.263678252006363	1/2.46695141192899	1/0.893774807719991	0.00297000000000000];
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

YTick_2D = 0:0.2:1;
XTick_2D = 0:1:6;
ZTick_2D = 0:300:1500;

YTickLabel_2D = {'0','0.2','0.4','0.6','0.8','1.0'};
XTickLabel_2D = {'0','1','2','3','4','5','6'};
ZTickLabel_2D = {'0','300','600','900','1200','1500'};

set(gca,'FontSize',10,'XTick',XTick_2D,'YTick',YTick_2D,'ZTick',ZTick_2D,'XTickLabel',XTickLabel_2D,'YTickLabel',YTickLabel_2D,'ZTickLabel',ZTickLabel_2D)
view(-70, 40); 

%% plot 2-params curve (Rc, C3)
load("contour T = T0-7 ~ T0+25, A1 = A10.mat")
curve = solution;

curve_params = zeros(size(curve,2),size(curve(1).params,2));

for i = 1:size(curve,2)
    curve_params(i,:) = curve(i).params;
end

figure
hold on

plot(1./curve_params(:,6),curve_params(:,1),LineWidth=lineWidth*2)

for i = 1:3
    para = params_plot{i};
    plot(para(6),para(1),MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
    plot(para(6),para(1),MarkerSize=mksize,Marker=".",Color=[0,0,0],LineStyle='--',LineWidth=lineWidth)
end

grid on
box on


axis([1 1.2 0 1500])

xlabel('{\itC}_3 (\muF)','Fontname','Arial')
ylabel('{\itR_C} (k\Omega)','Fontname','Arial')

XTick_2D = 1:0.05:1.2;
YTick_2D = 0:300:1500;

XTickLabel_2D = {'1.00','1,05','1.10','1.15','1.20'};
YTickLabel_2D = {'0','300','600','900','1200','1500'};

set(gca,'FontSize',10,'XTick',XTick_2D,'YTick',YTick_2D,'XTickLabel',XTickLabel_2D,'YTickLabel',YTickLabel_2D)
% view(-65, 40); 