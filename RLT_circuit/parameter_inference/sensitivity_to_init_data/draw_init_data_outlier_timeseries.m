clear
clc
%%
scaleLevels = [1.25];
repeatCount = 2;
indexPV = 1;
outputDirName = 'outlier_timeseries';
axisConfig.xLimit = [0 20];
axisConfig.xTicks = 0:5:20;
axisConfig.yLimit = [0 3];
axisConfig.yTicks = 0:0.6:3;
axisConfig.xLabel = 'Time (ms)';
axisConfig.yLabel = 'Voltage (V)';
axisConfig.fontSize = 10;

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
resultsRootDir = fullfile(scriptDir, 'results');
outputDir = fullfile(resultsRootDir, outputDirName);

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');
addpath(parameterInferenceDir, '-begin');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

learnedData = load(fullfile(circuitDir, 'learnedData_ODE.mat'), 'TS');
referenceData.t = learnedData.TS{1};
referenceData.y = learnedData.TS{2};
allOrbitL2 = [];

for scaleIdx = 1:numel(scaleLevels)
    statsData = load(fullfile(resultsRootDir, scale_dir_name(scaleLevels(scaleIdx)), ...
        'successful_params_stats.mat'), 'orbitRelativeL2', 'orbitSuccess');
    values = statsData.orbitRelativeL2(statsData.orbitSuccess);
    allOrbitL2 = [allOrbitL2; values(:)]; %#ok<AGROW>
end

[~, outlierThreshold] = largest_gap_cut(allOrbitL2, true);
fprintf('Outlier threshold from largest log-gap: %.6g\n', outlierThreshold);

outlierCount = 0;
for scaleIdx = 1:numel(scaleLevels)
    scale = scaleLevels(scaleIdx);
    scaleDirName = scale_dir_name(scale);
    resultScaleDir = fullfile(resultsRootDir, scaleDirName);
    statsData = load(fullfile(resultScaleDir, 'successful_params_stats.mat'));

    outlierIdx = find(statsData.orbitSuccess & statsData.orbitRelativeL2 >= outlierThreshold);
    for k = 1:numel(outlierIdx)
        sampleIdx = outlierIdx(k);
        outlierCount = outlierCount + 1;
        fileName = statsData.successFileNames{sampleIdx};
        finalOrbit = statsData.successFinalOrbit{sampleIdx};

        [fig, ax] = draw_timeseries_comparison( ...
            finalOrbit.t, finalOrbit.y, ...
            referenceData.t, referenceData.y, ...
            repeatCount, indexPV, ...
            scale, fileName, statsData.orbitRelativeL2(sampleIdx), axisConfig);

        outputFile = fullfile(outputDir, sprintf( ...
            'outlier_timeseries_scale_%s_%s.png', ...
            strrep(sprintf('%.12g', scale), '.', 'p'), erase(fileName, '.mat')));
        exportgraphics(ax, outputFile, 'Resolution', 300);

        fprintf('Saved outlier %d: scale=%.6g, file=%s, orbitL2=%.6g, output=%s\n', ...
            outlierCount, scale, fileName, statsData.orbitRelativeL2(sampleIdx), outputFile);
    end
end

fprintf('Saved %d outlier time-series figures to %s\n', outlierCount, outputDir);

function [fig, ax] = draw_timeseries_comparison(tInferred, yInferred, tReference, yReference, ...
    repeatCount, indexPV, scale, fileName, orbitL2, axisConfig)
lineWidth = 1.5;
colors_l = {'#CBBBC1', '#E4B7BC', '#F5E4C8'};
colors_s = {'#551F33', '#BD4146', '#ECC68C'};

[tPlotInferred, yPlotInferred] = build_periodic_curve( ...
    tInferred, yInferred, indexPV, repeatCount);
[tPlotReference, yPlotReference] = build_periodic_curve( ...
    tReference, yReference, indexPV, repeatCount);

fig = figure();
ax = axes(fig);
hold(ax, 'on')
grid(ax, 'on')
box(ax, 'on')

lineHandles = gobjects(1, size(yPlotInferred, 2));
referenceLineHandles = gobjects(1, size(yPlotReference, 2));
for stateIdx = 1:size(yPlotInferred, 2)
    lineHandles(stateIdx) = plot(ax, tPlotInferred, yPlotInferred(:, stateIdx), '-', ...
        'Color', colors_l{stateIdx}, 'LineWidth', lineWidth);
    referenceLineHandles(stateIdx) = plot(ax, ...
        tPlotReference, yPlotReference(:, stateIdx), '--', ...
        'Color', colors_s{stateIdx}, 'LineWidth', lineWidth);
end

axis(ax, [axisConfig.xLimit, axisConfig.yLimit]);

xlabel(ax, axisConfig.xLabel, 'FontName', 'Arial');
ylabel(ax, axisConfig.yLabel, 'FontName', 'Arial');
title(ax, sprintf('scale=%.3g, %s, orbit L2=%.4g', ...
    scale, erase(fileName, '.mat'), orbitL2), ...
    'FontName', 'Arial', 'FontWeight', 'normal', 'Interpreter', 'none');
set(ax, 'FontName', 'Arial', 'FontSize', axisConfig.fontSize, ...
    'XTick', axisConfig.xTicks, 'YTick', axisConfig.yTicks, ...
    'XTickLabel', compose('%.3g', axisConfig.xTicks), ...
    'YTickLabel', compose('%.3g', axisConfig.yTicks));
end

function [tPlot, yPlot] = build_periodic_curve(t, y, indexPV, repeatCount)
t = t(:);
[~, idxMax] = max(y(:, indexPV));

y = [y(idxMax:end-1, :); y(1:idxMax, :)];
t = [t(idxMax:end-1) - t(idxMax); t(1:idxMax) + t(end) - t(idxMax)];

tPlot = [];
yPlot = [];
tStart = 0;
for repeatIdx = 1:repeatCount
    if repeatIdx == 1
        tSeg = t + tStart;
        ySeg = y;
    else
        tSeg = t(2:end) + tStart;
        ySeg = y(2:end, :);
    end
    tPlot = [tPlot; tSeg]; %#ok<AGROW>
    yPlot = [yPlot; ySeg]; %#ok<AGROW>
    tStart = tPlot(end);
end
end

function [lowerCut, upperCut] = largest_gap_cut(values, useLogScale)
if useLogScale
    sortedValues = sort(log10(values(values > 0)));
    [~, gapIdx] = max(diff(sortedValues));
    lowerCut = 10 ^ sortedValues(gapIdx);
    upperCut = 10 ^ sortedValues(gapIdx + 1);
else
    sortedValues = sort(values);
    [~, gapIdx] = max(diff(sortedValues));
    lowerCut = sortedValues(gapIdx);
    upperCut = sortedValues(gapIdx + 1);
end
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end
