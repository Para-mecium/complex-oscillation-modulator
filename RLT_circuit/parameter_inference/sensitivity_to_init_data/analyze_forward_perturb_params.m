clear
clc

scriptDir = fileparts(mfilename('fullpath'));
summaryFile = fullfile(scriptDir, 'forward_perturbation_summary.mat');
saveFile = fullfile(scriptDir, 'forward_perturbation_stats.mat');

data = load(summaryFile);
baselineResult = data.baselineResult;
perturbationResults = data.perturbationResults;

numScenarios = numel(perturbationResults);
scenarioLabels = cell(numScenarios, 1);
periodRelChange = NaN(numScenarios, 1);
varAmpRelChange = NaN(numScenarios, 3);
forwardSuccess = false(numScenarios, 1);
perturbedPhysicalValue = NaN(numScenarios, 1);
baselinePhysicalValue = NaN(numScenarios, 1);

for i = 1:numScenarios
    directionSymbol = '+';
    if perturbationResults(i).direction < 0
        directionSymbol = '-';
    end

    scenarioLabels{i} = sprintf('%s %s%.0f%%', ...
        perturbationResults(i).paramName, directionSymbol, 100 * perturbationResults(i).level);
    periodRelChange(i) = perturbationResults(i).relChangePeriod;
    forwardSuccess(i) = perturbationResults(i).perturbedSuccess;
    perturbedPhysicalValue(i) = perturbationResults(i).perturbedPhysicalValue;
    baselinePhysicalValue(i) = perturbationResults(i).baselinePhysicalValue;

    if ~isempty(perturbationResults(i).relChangeVarAmp)
        varAmpRelChange(i, :) = reshape(perturbationResults(i).relChangeVarAmp, 1, []);
    end
end

statsSummary = struct();
statsSummary.numScenarios = numScenarios;
statsSummary.numSuccessfulForward = nnz(forwardSuccess);
statsSummary.forwardSuccessRate = nnz(forwardSuccess) / numScenarios;
statsSummary.meanAbsRelChangePeriod = mean(abs(periodRelChange(forwardSuccess)));
statsSummary.maxAbsRelChangePeriod = max(abs(periodRelChange(forwardSuccess)));
statsSummary.meanAbsRelChangeVarAmp = mean(abs(varAmpRelChange(forwardSuccess, :)), 1);
statsSummary.maxAbsRelChangeVarAmp = max(abs(varAmpRelChange(forwardSuccess, :)), [], 1);

save(saveFile, ...
    'baselineResult', 'scenarioLabels', 'forwardSuccess', ...
    'baselinePhysicalValue', 'perturbedPhysicalValue', ...
    'periodRelChange', 'varAmpRelChange', ...
    'statsSummary', 'perturbationResults');

fprintf('Saved forward perturbation statistics: %s\n', saveFile);
