lineWidth = 0.5;
mksize = 10;
stepsize = 10;
lambdaStepCap = 5e-3;

scriptDir = fileparts(mfilename('fullpath'));
pathFile = fullfile(scriptDir, 'params_modulation_path.mat');
targetFiles = { ...
    fullfile(scriptDir, 'period_target_1p0x.mat'), ...
    fullfile(scriptDir, 'period_target_1p5x.mat'), ...
    fullfile(scriptDir, 'period_target_2p0x.mat')};

addpath(scriptDir);
needGenerate = ~isfile(pathFile);
for i = 1:numel(targetFiles)
    needGenerate = needGenerate || ~isfile(targetFiles{i});
end
if needGenerate
    disp('Modulation path cache missing. Generate with Modulation_bi_target...');
    Modulation_bi_target('generatePathCache', true, 'lambdaStepCap', lambdaStepCap);
end

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

x3Bounds = padded_limits([c1_plot; marker_c1]);
y3Bounds = padded_limits([c2_plot; marker_c2]);
z3Bounds = padded_limits([rc_plot; marker_rc]);
x2Bounds = padded_limits([c3_plot; marker_c3]);
y2Bounds = padded_limits([rc_plot; marker_rc]);

%% plot 3-params curve (Rc, C1, C2)
figure
plot3(c1_plot, c2_plot, rc_plot, 'LineWidth', lineWidth * 2);
grid on
hold on
box on

surf([x3Bounds(1), x3Bounds(2); x3Bounds(1), x3Bounds(2)], ...
    marker_c2(1) * ones(2, 2), ...
    [z3Bounds(1), z3Bounds(1); z3Bounds(2), z3Bounds(2)], ...
    'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none')

for i = 1:numel(targetFiles)
    plot3([marker_c1(i), marker_c1(i)], [marker_c2(i), marker_c2(i)], ...
        [marker_rc(i), z3Bounds(1)], 'MarkerSize', mksize, 'Marker', '.', ...
        'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', lineWidth)
    plot3([marker_c1(i), marker_c1(i)], [marker_c2(i), y3Bounds(2)], ...
        [marker_rc(i), marker_rc(i)], 'MarkerSize', mksize, 'Marker', '.', ...
        'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', lineWidth)
end

axis([x3Bounds, y3Bounds, z3Bounds])

xlabel('{\itC}_1 (\muF)', 'Fontname', 'Arial')
ylabel('{\itC}_2 (\muF)', 'Fontname', 'Arial')
zlabel('{\itR_C} (k\Omega)', 'Fontname', 'Arial')
set(gca, 'FontSize', 10)
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

axis([x2Bounds, y2Bounds])

xlabel('{\itC}_3 (\muF)', 'Fontname', 'Arial')
ylabel('{\itR_C} (k\Omega)', 'Fontname', 'Arial')
set(gca, 'FontSize', 10)

%% local functions
function bounds = padded_limits(values)
values = double(values(:));
valueMin = min(values);
valueMax = max(values);
span = valueMax - valueMin;
scale = max(1, max(abs([valueMin, valueMax])));
pad = max(0.05 * span, 1e-3 * scale);
bounds = [valueMin - pad, valueMax + pad];
end
