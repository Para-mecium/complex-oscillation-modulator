clear
clc
%%
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
summaryFile = fullfile(scriptDir, 'params_inf_sensitivity_summary.mat');
saveEvery = 10;

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');

targetData = load(fullfile(parameterInferenceDir, 'initData_circuit.mat'));
targetPeriod = targetData.period;
targetVarAmp = reshape(targetData.varAmp, 1, []);

initDataFiles = dir(fullfile(scriptDir, 'init_data_files', 'initData_*.mat'));
if isempty(initDataFiles)
    error('params_inf_sensitivity_to_init_data:NoInitDataFiles', ...
        'No initData_*.mat files found in %s.', fullfile(scriptDir, 'init_data_files'));
end

[~, sortIdx] = sort({initDataFiles.name});
initDataFiles = initDataFiles(sortIdx);

runResults = repmat(struct( ...
    'fileName', '', ...
    'success', false, ...
    'finalLambda', NaN, ...
    'finalParams', [], ...
    'reason', '', ...
    'errorIdentifier', '', ...
    'errorMessage', ''), numel(initDataFiles), 1);

for k = 1:numel(initDataFiles)
    fileName = initDataFiles(k).name;
    filePath = fullfile(initDataFiles(k).folder, fileName);

    runResults(k).fileName = fileName;

    try
        initData = load(filePath);

        run(fullfile(circuitDir, 'System.m'));
        obs = [];

        params = initData.Parameters;
        t = initData.TS{1};
        TS_var = initData.TS{2};
        M = 50;
        derivatives = build_symbolic_derivatives(sys, obs, numel(params));

        PV.name = 'var';
        PV.idx = 1;
        State = state(obs, params, t, TS_var, M, PV);
        StateView = fmam_state_ops.solverViewFromState(State);

        items_per = struct;

        items_per(1).prop = 'p_Psi';
        items_per(1).idx = 1;
        items_per(1).target = targetPeriod / (2 * pi);

        items_per(2).prop = 'varAmp';
        items_per(2).idx = 1;
        items_per(2).target = targetVarAmp(1);

        items_per(3).prop = 'varAmp';
        items_per(3).idx = 2;
        items_per(3).target = targetVarAmp(2);

        items_per(4).prop = 'varAmp';
        items_per(4).idx = 3;
        items_per(4).target = targetVarAmp(3);

        items_controlled = [1 4 5 6];
        errBound = 1e-6;
        continuationOptions = struct('initialLambdaStep', 0.01, 'predictorMode', 'constant');

        Modtask = FMAM_ODE(sys, obs, StateView, items_per, items_controlled, [], errBound, ...
            'derivatives', derivatives, 'continuationOptions', continuationOptions);
        Modtask.psiUpdateMode = true;
        Modtask.refreshPsiModeReferences();
        Modtask.needLog = false;
        % Modtask.verbose = false;

        Modtask.fit();
        Modtask.step();
        Modtask.errBound = 1e-12;
        Modtask.fit();

        finalLambda = Modtask.continuationStatus.lambda;
        success = abs(finalLambda - 1) <= 1e-12;

        runResults(k).success = success;
        runResults(k).finalLambda = finalLambda;
        if success
            finalView = Modtask.exportSolverView();
            runResults(k).finalParams = reshape(finalView.params, 1, []);
        end
        runResults(k).reason = char(string(Modtask.continuationStatus.reason));

        fprintf('Processed %s: success=%d, lambda=%.12g\n', ...
            fileName, success, finalLambda);
    catch ME
        runResults(k).success = false;
        runResults(k).finalLambda = NaN;
        runResults(k).reason = 'exception';
        runResults(k).errorIdentifier = char(string(ME.identifier));
        runResults(k).errorMessage = char(string(ME.message));

        fprintf('Processed %s: success=0, lambda=NaN, error=%s\n', ...
            fileName, ME.identifier);
    end

    if mod(k, saveEvery) == 0
        numProcessed = k;
        partialRunResults = runResults(1:k);
        numSuccess = nnz([partialRunResults.success]);
        successRate = numSuccess / numProcessed;
        saveFile = fullfile(scriptDir, ...
            sprintf('params_inf_sensitivity_summary_%03d.mat', k));
        save(saveFile, 'numProcessed', 'numSuccess', 'successRate', 'partialRunResults');
        fprintf('Saved partial summary: %s\n', saveFile);
    end
end

numTotal = numel(runResults);
numSuccess = nnz([runResults.success]);
successRate = numSuccess / numTotal;

save(summaryFile, 'numTotal', 'numSuccess', 'successRate', 'runResults');

fprintf('Processed %d files, %d succeeded, success rate = %.6f\n', ...
    numTotal, numSuccess, successRate);
fprintf('Saved summary: %s\n', summaryFile);
