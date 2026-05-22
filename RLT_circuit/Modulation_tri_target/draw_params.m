clear
clc
%%
lineWidth = 0.5;
mksize = 10;
stepsize = 10;

scriptDir = fileparts(mfilename('fullpath'));
pathFile = fullfile(scriptDir, 'params_modulation_path.mat');
targetFiles = { ...
    fullfile(scriptDir, 'period_target_1p0x.mat'), ...
    fullfile(scriptDir, 'period_target_1p5x.mat'), ...
    fullfile(scriptDir, 'period_target_2p0x.mat')};

pathCache = load(pathFile);
curve_params = double(pathCache.curve_params);
params_start = reshape(double(pathCache.params_start), 1, []);
params_end = reshape(double(pathCache.params_end), 1, []);

if any(abs(curve_params(1, :) - params_start) > 1e-12)
    curve_params = [params_start; curve_params];
else
    curve_params(1, :) = params_start;
end
if any(abs(curve_params(end, :) - params_end) > 1e-12)
    curve_params = [curve_params; params_end];
end

sampleIdx = 1:stepsize:size(curve_params, 1);
if sampleIdx(end) ~= size(curve_params, 1)
    sampleIdx = [sampleIdx, size(curve_params, 1)];
end
curve_plot = curve_params(sampleIdx, :);

marker_params = zeros(numel(targetFiles), size(curve_params, 2));
for i = 1:numel(targetFiles)
    targetData = load(targetFiles{i});
    marker_params(i, :) = reshape(double(targetData.Parameters), 1, []);
end

rc_plot = curve_plot(:, 1);
c1_plot = 1 ./ curve_plot(:, 4);
c2_plot = 1 ./ curve_plot(:, 5);
c3_plot = 1 ./ curve_plot(:, 6);

marker_rc = marker_params(:, 1);
marker_c1 = 1 ./ marker_params(:, 4);
marker_c2 = 1 ./ marker_params(:, 5);
marker_c3 = 1 ./ marker_params(:, 6);

%% plot 3-params curve (Rc, C1, C2)
figure
plot3(c1_plot, c2_plot, rc_plot, 'LineWidth', lineWidth * 2);
grid on
hold on
box on

for i = 1:numel(targetFiles)
    plot3([marker_c1(i), marker_c1(i)], [marker_c2(i), marker_c2(i)], ...
        [marker_rc(i), -1], 'MarkerSize', mksize, 'Marker', '.', ...
        'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', lineWidth)
    plot3([marker_c1(i), marker_c1(i)], [marker_c2(i), 2], ...
        [marker_rc(i), marker_rc(i)], 'MarkerSize', mksize, 'Marker', '.', ...
        'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', lineWidth)
end

axis([1, 6, 0, 1, 0, 1500])

xlabel('{\itC}_1 (\muF)', 'Fontname', 'Arial')
ylabel('{\itC}_2 (\muF)', 'Fontname', 'Arial')
zlabel('{\itR_C} (k\Omega)', 'Fontname', 'Arial')

XTick_3D = 0:1:6;
YTick_3D = 0:0.2:1;
ZTick_3D = 0:300:1500;
XTickLabel_3D = {'0', '1', '2', '3', '4', '5', '6'};
YTickLabel_3D = {'0', '0.2', '0.4', '0.6', '0.8', '1.0'};
ZTickLabel_3D = {'0', '300', '600', '900', '1200', '1500'};

set(gca, 'FontSize', 10, 'XTick', XTick_3D, 'YTick', YTick_3D, ...
    'ZTick', ZTick_3D, 'XTickLabel', XTickLabel_3D, ...
    'YTickLabel', YTickLabel_3D, 'ZTickLabel', ZTickLabel_3D)
view(-70, 40);

%% plot 2-params curve (Rc, C3)
figure
plot(c3_plot, rc_plot, 'LineWidth', lineWidth * 2)
grid on
box on
hold on

for i = 1:numel(targetFiles)
    plot(marker_c3(i), marker_rc(i), 'MarkerSize', mksize, 'Marker', '.', ...
        'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', lineWidth)
end

axis([1, 1.2, 900, 1500])

xlabel('{\itC}_3 (\muF)', 'Fontname', 'Arial')
ylabel('{\itR_C} (k\Omega)', 'Fontname', 'Arial')

XTick_2D = 1:0.05:1.2;
YTick_2D = 0:300:1500;

XTickLabel_2D = {'1.00', '1.05', '1.10', '1.15', '1.20'};
YTickLabel_2D = {'0', '300', '600', '900', '1200', '1500'};

set(gca, 'FontSize', 10, 'XTick', XTick_2D, 'YTick', YTick_2D, ...
    'XTickLabel', XTickLabel_2D, 'YTickLabel', YTickLabel_2D)
