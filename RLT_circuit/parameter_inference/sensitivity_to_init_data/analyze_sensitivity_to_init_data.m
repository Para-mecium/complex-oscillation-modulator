clear
clc

scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
initDataRootDir = fullfile(scriptDir, 'init_data_files');
resultsRootDir = fullfile(scriptDir, 'results');
scaleLevels = [0.5 0.75 1 1.25 1.5];

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');
addpath(parameterInferenceDir, '-begin');

targetDataFile = fullfile(parameterInferenceDir, 'initData_circuit.mat');
targetData = load(targetDataFile);
targetOrbit = struct( ...
    't', targetData.t(:), ...
    'y', targetData.y, ...
    'period', targetData.period);
lossOptions = struct( ...
    'name', 'relative_l2_orbit', ...
    'compareNumPoints', size(targetData.y, 1), ...
    'phaseAlignment', true, ...
    'periodWeight', 0);

analysisByScale = repmat(struct( ...
    'scale', NaN, ...
    'resultScaleDir', '', ...
    'saveFile', '', ...
    'numTotal', NaN, ...
    'numSuccess', NaN, ...
    'successRate', NaN, ...
    'successFileNames', {{}}, ...
    'successSeeds', [], ...
    'successFinalParams', [], ...
    'successPhysicalParams', [], ...
    'successRuntimeFMAM', [], ...
    'successFinalOrbit', {{}}, ...
    'orbitRelativeL2', [], ...
    'orbitSuccess', [], ...
    'statsSummary', [], ...
    'runtimeSummaryFMAM', [], ...
    'orbitL2Summary', [], ...
    'corrMatrix', []), numel(scaleLevels), 1);
allOrbitRelativeL2 = zeros(0, 1);

for scaleIdx = 1:numel(scaleLevels)
    scale = scaleLevels(scaleIdx);
    initDataScaleDir = fullfile(initDataRootDir, scale_dir_name(scale));
    resultScaleDir = fullfile(resultsRootDir, scale_dir_name(scale));
    summaryFile = fullfile(resultScaleDir, 'params_inf_sensitivity_summary.mat');
    saveFile = fullfile(resultScaleDir, 'successful_params_stats.mat');

    summaryData = load(summaryFile);
    runResults = summaryData.runResults;
    numTotal = summaryData.numTotal;
    numSuccess = summaryData.numSuccess;
    successRate = summaryData.successRate;

    successIdx = find([runResults.success]);
    successResults = runResults(successIdx);
    successFileNames = {successResults.fileName}.';
    successFinalParams = vertcat(successResults.finalParams);
    successRuntimeFMAM = [successResults.runtimeFMAM].';
    successFinalOrbit = {successResults.finalOrbit}.';
    successSeeds = zeros(numel(successIdx), 1);
    orbitRelativeL2 = NaN(numel(successIdx), 1);
    orbitSuccess = [successResults.finalOrbitSuccess].';

    for k = 1:numel(successIdx)
        initData = load(fullfile(initDataScaleDir, successFileNames{k}));
        successSeeds(k) = initData.seed;
        if orbitSuccess(k)
            orbitRelativeL2(k) = loss_function( ...
                successFinalOrbit{k}, targetOrbit, lossOptions);
            orbitSuccess(k) = isfinite(orbitRelativeL2(k));
        end
    end

    successPhysicalParams = [ ...
        successFinalParams(:, 1), ...
        1 ./ successFinalParams(:, 4), ...
        1 ./ successFinalParams(:, 5), ...
        1 ./ successFinalParams(:, 6)];

    corrVarNames = {'R_C', 'C_1', 'C_2', 'C_3'};
    corrMatrix = corrcoef(successPhysicalParams);

    statsSummary = struct();
    statsSummary.R_C = summarize_column(successPhysicalParams(:, 1));
    statsSummary.C_1 = summarize_column(successPhysicalParams(:, 2));
    statsSummary.C_2 = summarize_column(successPhysicalParams(:, 3));
    statsSummary.C_3 = summarize_column(successPhysicalParams(:, 4));
    runtimeSummaryFMAM = summarize_column(successRuntimeFMAM);
    orbitL2Summary = summarize_column(orbitRelativeL2(orbitSuccess));

    analysisByScale(scaleIdx).scale = scale;
    analysisByScale(scaleIdx).resultScaleDir = resultScaleDir;
    analysisByScale(scaleIdx).saveFile = saveFile;
    analysisByScale(scaleIdx).numTotal = numTotal;
    analysisByScale(scaleIdx).numSuccess = numSuccess;
    analysisByScale(scaleIdx).successRate = successRate;
    analysisByScale(scaleIdx).successFileNames = successFileNames;
    analysisByScale(scaleIdx).successSeeds = successSeeds;
    analysisByScale(scaleIdx).successFinalParams = successFinalParams;
    analysisByScale(scaleIdx).successPhysicalParams = successPhysicalParams;
    analysisByScale(scaleIdx).successRuntimeFMAM = successRuntimeFMAM;
    analysisByScale(scaleIdx).successFinalOrbit = successFinalOrbit;
    analysisByScale(scaleIdx).orbitRelativeL2 = orbitRelativeL2;
    analysisByScale(scaleIdx).orbitSuccess = orbitSuccess;
    analysisByScale(scaleIdx).statsSummary = statsSummary;
    analysisByScale(scaleIdx).runtimeSummaryFMAM = runtimeSummaryFMAM;
    analysisByScale(scaleIdx).orbitL2Summary = orbitL2Summary;
    analysisByScale(scaleIdx).corrMatrix = corrMatrix;
    allOrbitRelativeL2 = [allOrbitRelativeL2; orbitRelativeL2(orbitSuccess)]; %#ok<AGROW>
end

[corrInlierCut, corrOutlierCut] = largest_gap_cut(allOrbitRelativeL2, true);

for scaleIdx = 1:numel(scaleLevels)
    scaleData = analysisByScale(scaleIdx);

    scale = scaleData.scale;
    saveFile = scaleData.saveFile;
    numTotal = scaleData.numTotal;
    numSuccess = scaleData.numSuccess;
    successRate = scaleData.successRate;
    successFileNames = scaleData.successFileNames;
    successSeeds = scaleData.successSeeds;
    successFinalParams = scaleData.successFinalParams;
    successPhysicalParams = scaleData.successPhysicalParams;
    successRuntimeFMAM = scaleData.successRuntimeFMAM;
    runtimeSummaryFMAM = scaleData.runtimeSummaryFMAM;
    successFinalOrbit = scaleData.successFinalOrbit;
    orbitRelativeL2 = scaleData.orbitRelativeL2;
    orbitSuccess = scaleData.orbitSuccess;
    orbitL2Summary = scaleData.orbitL2Summary;
    statsSummary = scaleData.statsSummary;
    corrMatrix = scaleData.corrMatrix;
    corrVarNames = {'R_C', 'C_1', 'C_2', 'C_3'};

    corrInlierMask = orbitSuccess & orbitRelativeL2 <= corrInlierCut;
    corrOutlierMask = orbitSuccess & orbitRelativeL2 >= corrOutlierCut;
    filteredSuccessPhysicalParams = successPhysicalParams(corrInlierMask, :);
    corrMatrixFiltered = corrcoef(filteredSuccessPhysicalParams);
    corrFilter = struct( ...
        'metric', 'orbit_l2', ...
        'useLogScale', true, ...
        'lowerCut', corrInlierCut, ...
        'upperCut', corrOutlierCut, ...
        'numMetricValues', nnz(orbitSuccess), ...
        'numInliers', nnz(corrInlierMask), ...
        'numOutliers', nnz(corrOutlierMask), ...
        'inlierMask', corrInlierMask, ...
        'outlierMask', corrOutlierMask);

    save(saveFile, ...
        'scale', 'targetDataFile', 'lossOptions', ...
        'numTotal', 'numSuccess', 'successRate', ...
        'successFileNames', 'successSeeds', ...
        'successFinalParams', 'successPhysicalParams', ...
        'successRuntimeFMAM', 'runtimeSummaryFMAM', ...
        'successFinalOrbit', 'orbitRelativeL2', 'orbitSuccess', 'orbitL2Summary', ...
        'statsSummary', 'corrMatrix', 'corrMatrixFiltered', ...
        'filteredSuccessPhysicalParams', 'corrFilter', 'corrVarNames');

    fprintf('scale=%.6g: saved successful parameter statistics: %s\n', ...
        scale, saveFile);
    fprintf(['scale=%.6g: successful samples: %d / %d (%.6f), ' ...
        'orbitSuccess=%d/%d, medianOrbitL2=%.6g, corrInliers=%d, corrOutliers=%d\n'], ...
        scale, numSuccess, numTotal, successRate, ...
        nnz(orbitSuccess), numel(orbitSuccess), orbitL2Summary.median, ...
        corrFilter.numInliers, corrFilter.numOutliers);
end

function summary = summarize_column(x)
summary = struct();
summary.mean = mean(x);
summary.median = median(x);
summary.std = std(x, 0);
summary.iqr = iqr(x);
summary.min = min(x);
summary.max = max(x);
summary.cv = summary.std / summary.mean;
summary.q25 = prctile(x, 25);
summary.q75 = prctile(x, 75);
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
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
