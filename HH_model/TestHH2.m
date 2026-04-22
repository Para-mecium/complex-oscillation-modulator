% 霍奇金-赫胥黎 (Hodgkin-Huxley) 模型：电流变化速率对双稳态系统的影响
% 左侧：电流极缓慢减小，系统绝热地跟随极限环，最终变为大振幅振荡（不消失）
% 右侧：电流在特定相位突然减小，系统状态落入平衡点吸引域，振荡消失

clear; clc; 

%% 1. 模型参数及可调电流定义
p.C = 1.0; 
p.gNa = 120.0; p.ENa = 50.0;
p.gK = 36.0;   p.EK = -77.0;
p.gL = 0.3;    p.EL = -54.387;
p.G = 0;

% ==========================================
% 在这里调节高电流和低电流的值
I_high = 120;  % 初始高电流（小振幅高频振荡）
I_low = 7;     % 最终低电流（具有双稳态：静息平衡点和大振幅极限环）
% ==========================================

y0 = [0, 0.5, 0.2, 0.5]; % 初始状态
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);

%% 2. 左侧情形：电流极缓慢减小 (Adiabatic change)
% 为了保证不脱离极限环，下降过程必须足够缓慢，这里设置总时间为 800ms
T_slow_total = 1000; 
% 50ms前保持 I_high，50ms到650ms缓慢降至 I_low，之后保持 I_low
I_slow = @(t) (t < 50) * I_high + ...
              (t >= 50 & t <= 650) .* (I_high - (I_high - I_low) * (t - 50) / 600) + ...
              (t > 650) * I_low;

fprintf('正在计算左侧缓慢变化情形 (需要较长时间以保证绝热性)...\n');
[t_slow, y_slow] = ode15s(@(t,y) hh_ode(t, y, I_slow(t), p), [0, T_slow_total], y0, options);
V_slow = y_slow(:, 1);
I_slow_val = I_slow(t_slow);

%% 3. 右侧情形：寻找合适的突变时刻 t_drop，使系统落入静息吸引域
T_step_total = 1000; % 右侧总模拟时间不需要太长

% 先在 I_high 下运行 500ms，让系统完全稳定在小极限环上
[t_init, y_init] = ode15s(@(t,y) hh_ode(t, y, I_high, p), [0, 500], y0, options);
y_steady = y_init(end, :);

% 在 500ms 到 505ms 之间（包含多个高频振荡周期）进行高精度扫描
t_search = linspace(500, 505, 300);
[~, y_search] = ode15s(@(t,y) hh_ode(t, y, I_high, p), t_search, y_steady, options);

t_drop = 500; % 初始化
found = false;
for i = 1:length(t_search)
    y_test_start = y_search(i, :);
    % 从该点开始，以 I_low 运行 50ms 进行测试
    [~, y_test] = ode15s(@(t,y) hh_ode(t, y, I_low, p), [0, 50], y_test_start, options);
    
    % 检查是否产生动作电位：如果后半段的电压始终低于 -20mV，说明没有放电，成功落入吸引域
    if max(y_test(20:end, 1)) < -20 
        t_drop = t_search(i);
        found = true;
        break;
    end
end

if found
    fprintf('找到合适的突降时刻: t_drop = %.3f ms\n', t_drop);
else
    fprintf('警告: 未找到完美的突降时刻，使用默认值。\n');
end

% 定义阶跃变化的电流函数
I_step = @(t) (t < t_drop) * I_high + (t >= t_drop) * I_low;

[t_step, y_step] = ode15s(@(t,y) hh_ode(t, y, I_step(t), p), [0, T_step_total], y0, options);
V_step = y_step(:, 1);
I_step_val = I_step(t_step);

%% 4. 可视化
figure('Position', [50, 100, 1200, 600], 'Name', 'Bistability and Parameter Change Rate');

% 电流图的 Y 轴范围自适应
I_ylim = [min(0, I_low - 10), I_high + 10];

% --- 左侧：缓慢减小 ---
subplot(2, 2, 1);
plot(t_slow, V_slow, 'b-', 'LineWidth', 1.2);
title('Neuron 1: Slow Current Decrease (Stays on Limit Cycle)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (mV)', 'FontSize', 11);
ylim([-90, 50]); grid on;
xlim([0, T_slow_total]);

subplot(2, 2, 3);
plot(t_slow, I_slow_val, 'r-', 'LineWidth', 2);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Current (\muA/cm^2)', 'FontSize', 11);
ylim(I_ylim); grid on;
xlim([0, T_slow_total]);

% --- 右侧：突然减小 ---
subplot(2, 2, 2);
plot(t_step, V_step, 'b-', 'LineWidth', 1.2);
title('Neuron 2: Sudden Current Drop (Falls into Resting Basin)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (mV)', 'FontSize', 11);
ylim([-90, 50]); grid on;
xlim([0, T_step_total]);
% 标注突变时刻
xline(t_drop, 'k--', sprintf('Drop to I=%g', I_low), 'LabelVerticalAlignment', 'bottom', 'FontSize', 10);

subplot(2, 2, 4);
plot(t_step, I_step_val, 'r-', 'LineWidth', 2);
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Current (\muA/cm^2)', 'FontSize', 11);
ylim(I_ylim); grid on;
xlim([0, T_step_total]);

sgtitle('Sensitivity to Parameter Change Rate in Hodgkin-Huxley Model', 'FontSize', 16, 'FontWeight', 'bold');

%% ================= 辅助函数 ================= %%

% HH 模型的常微分方程
function dydt = hh_ode(~, y, I, p)
    V = y(1); m = y(2); h = y(3); n = y(4);
    
    [am, bm, ah, bh, an, bn] = get_rates(V);
    
    dVdt = (I - p.gNa*m^3*h*(V - p.ENa) - p.gK*n^4*(V - p.EK) - p.gL*(V - p.EL)) / p.C;
    dmdt = am*(1-m) - bm*m;
    dhdt = ah*(1-h) - bh*h;
    dndt = an*(1-n) - bn*n;
    
    dydt = [dVdt; dmdt; dhdt; dndt];
end

% 计算速率常数
function [am, bm, ah, bh, an, bn] = get_rates(V)
    am = 0.1 * (V + 40) / (1 - exp(-(V + 40)/10));
    bm = 4 * exp(-(V + 65)/18);
    ah = 0.07 * exp(-(V + 65)/20);
    bh = 1 / (1 + exp(-(V + 35)/10));
    an = 0.01 * (V + 55) / (1 - exp(-(V + 55)/10));
    bn = 0.125 * exp(-(V + 65)/80);
    
    if isnan(am), am = 1; end
    if isnan(an), an = 0.1; end
end