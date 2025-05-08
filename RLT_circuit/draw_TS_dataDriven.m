clear
clc
%% pre-setting
linewidth = 1;
index_PV = 1;
mksize = 25;
coe_deeper = 0.8;
% colors_l = {'#F7C2CD','#A6DAEF','#B0D9A5'};
% colors_s = {'#E26472','#6270B7','#077535'};

colors_l = {'#CBBBC1','#E4B7BC','#F5E4C8'};
colors_s = {'#551F33','#BD4146','#ECC68C'};
% 
% colors_l = {[191 203 249]/255,[248 224 189]/255,[224 243 220]/255};
% colors_s = {[156 168 228]/255,[235 194 140]/255,[173 214 165]/255};


%% ODE plot
cls = 'init';
t_plot = [];
y_plot = [];
tStart = 0;

if strcmpi(cls,'init')
    load("initData_ODE.mat")
elseif strcmpi(cls,'learned')
    load("learnedData_ODE.mat")
end

t = TS{1};
y = TS{2};

[~,idx_max] = max(y(:,index_PV));
y = [y(idx_max:end-1,:);y(1:idx_max,:)]/0.06;
t = [t(idx_max:end-1)-t(idx_max);t(1:idx_max)+t(end)-t(idx_max)];
if strcmpi(cls,'init')
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    y_plot = [y_plot;y;y];
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    y_plot = [y_plot;y;y];
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    y_plot = [y_plot;y;y];
elseif strcmpi(cls,'learned')
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    t_plot = [t_plot;t+tStart];
    tStart = t_plot(end);
    y_plot = [y_plot;y;y];
end

figure
hold on

plt_l = plot(t_plot,y_plot,LineStyle='-',LineWidth=linewidth);

grid on
box on

%% circuit plot
ptNum = 200;
stepsize = 4;
t_plot = [];
y_plot = [];
tStart = 0;

load("initData_circuit.mat")

[~,idx_max] = max(y(:,index_PV));
y = [y(idx_max:end-1,:);y(1:idx_max,:)]/0.06;
t = [t(idx_max:end-1)-t(idx_max);t(1:idx_max)+t(end)-t(idx_max)];

t_plot = [t_plot;t+tStart];
tStart = t_plot(end);
t_plot = [t_plot;t(2:end)+tStart];
tStart = t_plot(end);
y_plot = [y_plot;y;y(2:end,:)];

tend = t_plot(end);
tt = (0:1:ptNum-1)'/ptNum*tend;
for i = 1:3
    yy(:,i) = spline(t_plot,y_plot(:,i),tt);
end
t_plot = tt;
y_plot = yy;

plt_s = scatter(t_plot(1:stepsize:end),y_plot(1:stepsize:end,:),mksize,Marker="x");


xlabel('Time (a.u.)','FontName','Arial')
ylabel('Concentration (a.u.)','FontName','Arial')

XTick_2DTS = [0 5 10 15 20];
YTick_2DTS = [0 10 20 30 40 50];
XTickLabel_2DTS = {'0','5','10','15','20'};
YTickLabel_2DTS = {'0','10','20','30','40','50'};

set(gca,'FontSize',10,'XTick',XTick_2DTS,'YTick',YTick_2DTS,'XTickLabel',XTickLabel_2DTS,'YTickLabel',YTickLabel_2DTS)

for i = 1:3
    set(plt_l(i),'Color',colors_l{i})
    set(plt_s(i),'MarkerFaceColor',colors_s{i},'MarkerEdgeColor',colors_s{i})
end
axis([0 20 0 50])