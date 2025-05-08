function [t_period, y_period] = extractPeriod(t, y, method, varargin)
% EXTRACTCONSISTENTPERIOD 提取首尾一致且以最大值为起点的周期序列
%   [t_period, y_period] = extractConsistentPeriod(t, y) 自动检测周期并提取
%   [t_period, y_period] = extractConsistentPeriod(t, y, 'fixed', T) 使用已知周期T
%
% 特点:
%   1. 保证提取的周期序列首尾值一致(相位连续)
%   2. 以周期内最大值为起始点
%   3. 自动对齐多个通道(使用第一个通道确定相位)

% 参数检查
if nargin < 3
    method = 'auto';
end

% 固定周期模式
if strcmpi(method, 'fixed')
    if nargin < 4
        error('使用fixed模式需要提供周期T');
    end
    T = varargin{1};
else
    % 自动检测周期 ---------------------------------------------------
    % 使用第一个变量检测周期
    sig = y(:,1);
    
    % 方法1: 自相关分析
    [acf, lags] = xcorr(sig-mean(sig), 'coeff');
    lags = lags(lags>0);
    acf = acf(lags>0);
    
    % 寻找显著峰值
    [peaks, locs] = findpeaks(acf, 'MinPeakHeight', 0.5, 'MinPeakDistance', 10);
    
    if isempty(locs)
        error('无法检测到显著周期，请尝试手动指定周期');
    end
    
    % 计算峰值间距作为周期估计
    peak_intervals = diff(lags(locs));
    T = round(median(peak_intervals));
    
    % 方法2: 傅里叶变换验证
    Fs = 1/mean(diff(t));
    n = length(sig);
    Y = fft(sig);
    P2 = abs(Y/n);
    P1 = P2(1:floor(n/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:floor(n/2))/n;
    
    [~, idx] = max(P1(2:floor(n/2))); % 忽略直流分量
    T_fft = round(1/f(idx+1)*Fs);
    
    % 选择更可靠的估计
    if abs(T - T_fft)/min(T, T_fft) > 0.2
        warning('两种方法检测结果不一致，使用自相关结果');
    else
        T = round((T + T_fft)/2);
    end
end

% 相位对齐: 找到周期内最大值位置 -----------------------------------
% 使用第一个变量确定相位
[~, max_idx] = max(y(:,1));

% 确保不超过数据范围
start_idx = max_idx; % 转换为第一个周期内的位置


% 提取初始周期
t_candidate = t(start_idx:start_idx+T-1);
y_candidate = y(start_idx:start_idx+T-1, :);

% 首尾一致性检查与调整 -------------------------------------------
tol = 1e-2 * std(y(:,1)); % 允许的差异阈值

% 检查第一个变量的首尾一致性
while abs(y_candidate(1,1) - y_candidate(end,1)) > tol
    % 微调起点使首尾更接近
    if y_candidate(1,1) > y_candidate(end,1)
        start_idx = start_idx - 1;
    else
        start_idx = start_idx + 1;
    end
    
    % 边界检查
    if start_idx <= 0 || start_idx + T - 1 > length(t)
        warning('无法找到完美匹配的周期，使用最佳近似');
        break;
    end
    
    % 重新提取
    t_candidate = t(start_idx:start_idx+T-1);
    y_candidate = y(start_idx:start_idx+T-1, :);
end

% 最终提取
t_period = t_candidate;
y_period = y_candidate;

% 可视化验证 -----------------------------------------------------
figure;

% 原始信号与提取的周期
subplot(3,1,1);
plot(t, y(:,1)), hold on;
plot(t_period, y_period(:,1), 'linewidth', 2);
plot(t_period([1 end]), y_period([1 end],1), 'ro', 'MarkerSize', 8);
title('原始信号与提取周期 (红点显示首尾)');
legend('原始信号', '提取周期');

% 周期序列
subplot(3,1,2);
plot(t_period - t_period(1), y_period);
title('对齐后的周期序列');
xlabel('时间 (相对)');
ylabel('幅值');

% 首尾连接检查
subplot(3,1,3);
plot([y_period(:,1); y_period(1,1)], 'o-');
title('首尾连接检查');
xlabel('样本点');
ylabel('幅值');
end