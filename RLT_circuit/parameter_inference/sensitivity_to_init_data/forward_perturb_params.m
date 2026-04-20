clear
clc
%%
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
statsFile = fullfile(scriptDir, 'successful_params_stats.mat');
saveFile = fullfile(scriptDir, 'forward_perturbation_summary.mat');

perturbationLevels = [0.01, 0.05];
directionValues = [-1, 1];
physicalParamNames = {'R_C', 'C_1', 'C_2', 'C_3'};

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');

statsData = load(statsFile);
if isempty(statsData.successFinalParams)
    error('forward_perturb_params:NoSuccessfulParams', ...
        'No successful parameters found in %s.', statsFile);
end

initDataODE = load(fullfile(parameterInferenceDir, 'initData_ODE.mat'));
y0 = initDataODE.TS{2}(1, :).';

baselineParams = median(statsData.successFinalParams, 1);
baselinePhysicalParams = median(statsData.successPhysicalParams, 1);
baselineParams(1) = baselinePhysicalParams(1);
baselineParams(4) = 1 / baselinePhysicalParams(2);
baselineParams(5) = 1 / baselinePhysicalParams(3);
baselineParams(6) = 1 / baselinePhysicalParams(4);

numScenarios = numel(physicalParamNames) * numel(perturbationLevels) * numel(directionValues);

fprintf('Forward perturbation starts: median baseline, %d scenarios.\n', numScenarios);

baselineResult = struct( ...
    'baselineSuccess', false, ...
    'message', '', ...
    'baselineParams', baselineParams, ...
    'baselinePhysicalParams', baselinePhysicalParams, ...
    'baselinePeriod', NaN, ...
    'baselineVarAmp', []);

baselineOrbit = circuit_forward_orbit(baselineParams, y0, struct());
if ~baselineOrbit.success
    baselineResult.message = baselineOrbit.msg;
    fprintf('Baseline failed: %s\n', baselineOrbit.msg);
else
    baselineResult.baselineSuccess = true;
    baselineResult.baselinePeriod = baselineOrbit.features.period;
    baselineResult.baselineVarAmp = reshape(baselineOrbit.features.state.amplitude, 1, []);
    fprintf('Baseline success: period=%.6f\n', baselineResult.baselinePeriod);
end

baselineSuccess = baselineResult.baselineSuccess;
baselineMessage = baselineResult.message;
baselinePeriod = baselineResult.baselinePeriod;
baselineVarAmp = baselineResult.baselineVarAmp;

perturbationResults = repmat(struct( ...
    'paramName', '', ...
    'paramIndex', NaN, ...
    'level', NaN, ...
    'direction', NaN, ...
    'baselineSuccess', baselineSuccess, ...
    'perturbedSuccess', false, ...
    'message', '', ...
    'baselinePhysicalValue', NaN, ...
    'perturbedPhysicalValue', NaN, ...
    'baselinePeriod', NaN, ...
    'perturbedPeriod', NaN, ...
    'relChangePeriod', NaN, ...
    'baselineVarAmp', [], ...
    'perturbedVarAmp', [], ...
    'relChangeVarAmp', []), numScenarios, 1);

scenarioIdx = 0;
for paramIndex = 1:numel(physicalParamNames)
    for level = perturbationLevels
        for direction = directionValues
            scenarioIdx = scenarioIdx + 1;

            scenarioResult = perturbationResults(scenarioIdx);
            scenarioResult.paramName = physicalParamNames{paramIndex};
            scenarioResult.paramIndex = paramIndex;
            scenarioResult.level = level;
            scenarioResult.direction = direction;

            if ~baselineSuccess
                scenarioResult.message = baselineMessage;
                perturbationResults(scenarioIdx) = scenarioResult;
                continue
            end

            perturbedPhysicalParams = baselinePhysicalParams;
            perturbedPhysicalParams(paramIndex) = baselinePhysicalParams(paramIndex) * (1 + direction * level);

            perturbedParams = baselineParams;
            perturbedParams(1) = perturbedPhysicalParams(1);
            perturbedParams(4) = 1 / perturbedPhysicalParams(2);
            perturbedParams(5) = 1 / perturbedPhysicalParams(3);
            perturbedParams(6) = 1 / perturbedPhysicalParams(4);

            scenarioResult.baselinePhysicalValue = baselinePhysicalParams(paramIndex);
            scenarioResult.perturbedPhysicalValue = perturbedPhysicalParams(paramIndex);
            scenarioResult.baselinePeriod = baselinePeriod;
            scenarioResult.baselineVarAmp = baselineVarAmp;

            perturbedOrbit = circuit_forward_orbit(perturbedParams, y0, struct());
            if ~perturbedOrbit.success
                scenarioResult.message = perturbedOrbit.msg;
                perturbationResults(scenarioIdx) = scenarioResult;
                continue
            end

            perturbedVarAmp = reshape(perturbedOrbit.features.state.amplitude, 1, []);
            scenarioResult.perturbedSuccess = true;
            scenarioResult.perturbedPeriod = perturbedOrbit.features.period;
            scenarioResult.relChangePeriod = ...
                (perturbedOrbit.features.period - baselinePeriod) / baselinePeriod;
            scenarioResult.perturbedVarAmp = perturbedVarAmp;
            scenarioResult.relChangeVarAmp = ...
                (perturbedVarAmp - baselineVarAmp) ./ baselineVarAmp;
            perturbationResults(scenarioIdx) = scenarioResult;
        end
    end
end

numForwardSuccess = nnz([perturbationResults.perturbedSuccess]);
fprintf('Perturbation done: %d/%d forward runs succeeded.\n', ...
    numForwardSuccess, numScenarios);

groupSummary = repmat(struct( ...
    'paramName', '', ...
    'level', NaN, ...
    'direction', NaN, ...
    'numCases', 0, ...
    'numSuccessfulForward', 0, ...
    'forwardSuccessRate', NaN, ...
    'meanAbsRelChangePeriod', NaN, ...
    'maxAbsRelChangePeriod', NaN, ...
    'meanAbsRelChangeVarAmp', [], ...
    'maxAbsRelChangeVarAmp', []), numScenarios, 1);

for scenarioIdx = 1:numScenarios
    successMask = perturbationResults(scenarioIdx).perturbedSuccess;

    groupSummary(scenarioIdx).paramName = perturbationResults(scenarioIdx).paramName;
    groupSummary(scenarioIdx).level = perturbationResults(scenarioIdx).level;
    groupSummary(scenarioIdx).direction = perturbationResults(scenarioIdx).direction;
    groupSummary(scenarioIdx).numCases = 1;
    groupSummary(scenarioIdx).numSuccessfulForward = double(successMask);
    groupSummary(scenarioIdx).forwardSuccessRate = double(successMask);

    if successMask
        relPeriod = abs(perturbationResults(scenarioIdx).relChangePeriod);
        relVarAmp = abs(perturbationResults(scenarioIdx).relChangeVarAmp);
        groupSummary(scenarioIdx).meanAbsRelChangePeriod = relPeriod;
        groupSummary(scenarioIdx).maxAbsRelChangePeriod = relPeriod;
        groupSummary(scenarioIdx).meanAbsRelChangeVarAmp = relVarAmp;
        groupSummary(scenarioIdx).maxAbsRelChangeVarAmp = relVarAmp;
    end
end

save(saveFile, 'perturbationLevels', 'directionValues', 'physicalParamNames', ...
    'baselineResult', 'perturbationResults', 'groupSummary');

fprintf('Saved forward perturbation summary: %s\n', saveFile);
