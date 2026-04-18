lineWidth = 0.5;
mksize = 10;
stepsize = 1;

scriptDir = fileparts(mfilename('fullpath'));
plotDataFile = fullfile(scriptDir, "params_modulation_path.mat");

dataCache = load(plotDataFile);
curve_params = dataCache.curve_params;
params_start = dataCache.params_start;
params_end = dataCache.params_end;

params_start = reshape(params_start, 1, []);
params_end = reshape(params_end, 1, []);

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

rc_all = curve_params(:, 1);
c1_all = 1 ./ curve_params(:, 4);
c2_all = 1 ./ curve_params(:, 5);
c3_all = 1 ./ curve_params(:, 6);

rc_plot = curve_plot(:, 1);
c1_plot = 1 ./ curve_plot(:, 4);
c2_plot = 1 ./ curve_plot(:, 5);
c3_plot = 1 ./ curve_plot(:, 6);

start3 = [1 / params_start(4), 1 / params_start(5), params_start(1)];
end3 = [1 / params_end(4), 1 / params_end(5), params_end(1)];

%% plot 3-params curve (Rc, C1, C2)
figure
plot3(c1_plot, c2_plot, rc_plot, 'LineWidth', lineWidth * 2);
grid on
hold on
box on

plot3(start3(1), start3(2), start3(3), 'k.', 'MarkerSize', mksize);
plot3(end3(1), end3(2), end3(3), 'k.', 'MarkerSize', mksize);

xPad = max(1e-6, 0.05 * (max(c1_all) - min(c1_all)));
yPad = max(1e-6, 0.05 * (max(c2_all) - min(c2_all)));
zPad = max(1e-6, 0.05 * (max(rc_all) - min(rc_all)));
axis([min(c1_all) - xPad, max(c1_all) + xPad, ...
      min(c2_all) - yPad, max(c2_all) + yPad, ...
      min(rc_all) - zPad, max(rc_all) + zPad]);

xlabel('{\itC}_1 (\muF)', 'Fontname', 'Arial')
ylabel('{\itC}_2 (\muF)', 'Fontname', 'Arial')
zlabel('{\itR_C} (k\Omega)', 'Fontname', 'Arial')
set(gca, 'FontSize', 10)
view(-70, 40);

%% plot 2-params curve (Rc, C3)
figure
plot(c3_plot, rc_plot, 'LineWidth', lineWidth * 2);
grid on
box on
hold on

plot(1 / params_start(6), params_start(1), 'k.', 'MarkerSize', mksize);
plot(1 / params_end(6), params_end(1), 'k.', 'MarkerSize', mksize);

xPad = max(1e-6, 0.05 * (max(c3_all) - min(c3_all)));
yPad = max(1e-6, 0.05 * (max(rc_all) - min(rc_all)));
axis([min(c3_all) - xPad, max(c3_all) + xPad, ...
      min(rc_all) - yPad, max(rc_all) + yPad]);

xlabel('{\itC}_3 (\muF)', 'Fontname', 'Arial')
ylabel('{\itR_C} (k\Omega)', 'Fontname', 'Arial')
set(gca, 'FontSize', 10)
