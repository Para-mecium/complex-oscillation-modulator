clear
clc

scriptDir = fileparts(mfilename('fullpath'));
summaryFile = fullfile(scriptDir, 'params_inf_sensitivity_summary.mat');
saveFile = fullfile(scriptDir, 'successful_params_stats.mat');

summaryData = load(summaryFile);
runResults = summaryData.runResults;
numTotal = summaryData.numTotal;
numSuccess = summaryData.numSuccess;
successRate = summaryData.successRate;

successIdx = find([runResults.success]);
if isempty(successIdx)
    error('analyze_successful_params:NoSuccessfulSamples', ...
        'No successful samples found in %s.', summaryFile);
end

successResults = runResults(successIdx);
successFileNames = {successResults.fileName}.';
successFinalParams = vertcat(successResults.finalParams);
successSeeds = zeros(numel(successIdx), 1);

for k = 1:numel(successIdx)
    initData = load(fullfile(scriptDir, 'init_data_files', successFileNames{k}));
    successSeeds(k) = initData.seed;
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

save(saveFile, ...
    'numTotal', 'numSuccess', 'successRate', ...
    'successFileNames', 'successSeeds', ...
    'successFinalParams', 'successPhysicalParams', ...
    'statsSummary', 'corrMatrix', 'corrVarNames');

fprintf('Saved successful parameter statistics: %s\n', saveFile);
fprintf('Successful samples: %d / %d (%.6f)\n', numSuccess, numTotal, successRate);

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
