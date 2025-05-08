clear
clc
%% pre-setting
index_PV = 1;
period_multiplier = 2; % 1, 1.5, 2

period_num = 6;
linewidth = 1;
ptNumAll = 90; 
ptNum = ptNumAll/period_num*period_multiplier;
mksize = 25;

startidx = 1;

colors_l = {'#CBBBC1','#E4B7BC','#F5E4C8'};
colors_s = {'#551F33','#BD4146','#ECC68C'};
%% ODE plot 
t_plot = [];
y_plot = [];
timeStart = 0;

load(['TS T = T0*' num2str(period_multiplier) ', A1 = A10.mat'])
t = TS{1};
y = TS{2}/0.06;


% [~,idx_min] = min(y(:,index_PV));
% y = [y(idx_min:end-1,:);y(1:idx_min,:)];
% t = [t(idx_min:end-1)-t(idx_min);t(1:idx_min)+t(end)-t(idx_min)];

count = 0;
while count < period_num
    count = count + 1;
    t_plot = [t_plot;t+timeStart];
    timeStart = t_plot(end);
    y_plot = [y_plot;y];
end

figure
plt_l = plot(t_plot,y_plot(:,index_PV),LineStyle='-',LineWidth=linewidth);
grid on
box on
hold on

%% circuit plot
TS_period_10 = readtable("1x_10.txt");
TS_period_10 = renamevars(TS_period_10,"x_____","t");
 
TS_period_15 = readtable("1.5x_10.txt");
TS_period_15 = renamevars(TS_period_15,"x_____","t");

TS_period_20 = readtable("2x_10.txt");
TS_period_20 = renamevars(TS_period_20,"x_____","t");

tabs = struct;
tabs.tab10 = TS_period_10;
tabs.tab15 = TS_period_15;
tabs.tab20 = TS_period_20;

tab = tabs.(['tab' num2str(10*period_multiplier)]);

t = tab.t(startidx:end)*1e3;
V1 = tab.VF1(startidx:end)/0.06;
V2 = tab.VF2(startidx:end)/0.06;
V3 = tab.VF3(startidx:end)/0.06;

y = [V1 V2 V3];

PointNum = size(t,1);
localmax = islocalmax(y(:,1));

index_max = [];
t_plot = [];
y_plot = [];

timeStart = 0;

for i = 1:PointNum
    if localmax(i) == 1
        index_max = [index_max i];
    end
end

index_start = index_max(end-1);
index_end = index_max(end);

t = t(index_start:index_end) - t(index_start);
y = y(index_start:index_end,:);

tend = t(end);
tstart = t(1);
tt = (0:1:ptNum)'/ptNum*(tend-tstart) + tstart;
for i = 1:3
    yy(:,i) = spline(t,y(:,i),tt);
end
t = tt;
y = yy;

count = 0;
while count < period_num
    count = count + 1;
    t_plot = [t_plot;t+timeStart];
    timeStart = t_plot(end);
    y_plot = [y_plot;y];
end

plt_s = scatter(t_plot,y_plot(:,index_PV),mksize,marker='x');

%%
T_0 = 10.6519;
xlabel('Time (a.u.)','FontName','Arial')
ylabel('Concentration (a.u.)','FontName','Arial')

count = 0;
XTick_2DTS = 0;
XTickLabel_2DTS = {'0'};
while count < period_num
    count = count + 1;
    XTick_2DTS(count) = count*T_0;
    if mod(count,2) == 0
        XTickLabel_2DTS{count} = [num2str(count) 'T_0'];
    else 
        XTickLabel_2DTS{count} = '';
    end
end

if index_PV == 1
    YTick_2DTS = [0 5 10 15 20 25];
    YTickLabel_2DTS = {'0','5','10','15','20','25'};
elseif index_PV == 2
    YTick_2DTS = [0 15 30 45 60 75];
    YTickLabel_2DTS = {'0','15','30','45','60','75'};
elseif index_PV == 3
    YTick_2DTS = [0 15 30 45 60 75];
    YTickLabel_2DTS = {'0','15','30','45','60','75'};
end
axis([0 T_0*period_num 0 YTick_2DTS(end)])

set(gca,'FontSize',10,'XTick',XTick_2DTS,'YTick',YTick_2DTS,'XTickLabel',XTickLabel_2DTS,'YTickLabel',YTickLabel_2DTS)
set(plt_l,'Color',colors_l{index_PV})
set(plt_s,'MarkerFaceColor',colors_s{index_PV},'MarkerEdgeColor',colors_s{index_PV})
