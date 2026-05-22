clear
clc
%%
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
initDataRootDir = fullfile(scriptDir, 'init_data_files');
resultsRootDir = fullfile(scriptDir, 'results');
scaleLevels = [0.5 0.75 1 1.25 1.5];
saveEvery = 10;

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');

targetData = load(fullfile(parameterInferenceDir, 'initData_circuit.mat'));
targetPeriod = targetData.period;
targetVarAmp = reshape(targetData.varAmp, 1, []);

for scaleIdx = 1:numel(scaleLevels)
    scale = scaleLevels(scaleIdx);
    initDataScaleDir = fullfile(initDataRootDir, scale_dir_name(scale));
    resultScaleDir = fullfile(resultsRootDir, scale_dir_name(scale));
    summaryFile = fullfile(resultScaleDir, 'params_inf_sensitivity_summary.mat');

    if ~isfolder(resultScaleDir)
        mkdir(resultScaleDir);
    end

    initDataFiles = dir(fullfile(initDataScaleDir, 'initData_*.mat'));
    [~, sortIdx] = sort({initDataFiles.name});
    initDataFiles = initDataFiles(sortIdx);

    runResults = repmat(struct( ...
        'fileName', '', ...
        'success', false, ...
        'finalLambda', NaN, ...
        'finalParams', [], ...
        'runtimeFMAM', NaN, ...
        'finalOrbit', [], ...
        'finalOrbitSuccess', false, ...
        'finalOrbitMessage', '', ...
        'reason', '', ...
        'errorIdentifier', '', ...
        'errorMessage', ''), numel(initDataFiles), 1);

    fprintf('Processing scale=%.6g, files=%d\n', scale, numel(initDataFiles));

    for k = 1:numel(initDataFiles)
        fileName = initDataFiles(k).name;
        filePath = fullfile(initDataFiles(k).folder, fileName);

        runResults(k).fileName = fileName;
        solveTimer = [];

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
            Modtask.verbose = false;

            solveTimer = tic;
            Modtask.fit();
            Modtask.step();
            Modtask.errBound = 1e-12;
            Modtask.fit();
            runResults(k).runtimeFMAM = toc(solveTimer);

            finalLambda = Modtask.continuationStatus.lambda;
            success = abs(finalLambda - 1) <= 1e-12;

            runResults(k).success = success;
            runResults(k).finalLambda = finalLambda;
            if success
                finalView = Modtask.exportSolverView();
                finalParams = reshape(finalView.params, 1, []);
                initialState = initData.TS{2}(1, :).';
                forwardResult = circuit_forward_orbit(finalParams, initialState, struct());

                runResults(k).finalParams = finalParams;
                runResults(k).finalOrbitSuccess = forwardResult.success;
                runResults(k).finalOrbitMessage = char(string(forwardResult.msg));
                if forwardResult.success
                    runResults(k).finalOrbit = forwardResult.orbit;
                end
            end
            runResults(k).reason = char(string(Modtask.continuationStatus.reason));

            fprintf(['Processed scale=%.6g %s: success=%d, lambda=%.12g, ' ...
                'runtimeFMAM=%.6g, orbitSuccess=%d\n'], ...
                scale, fileName, success, finalLambda, ...
                runResults(k).runtimeFMAM, runResults(k).finalOrbitSuccess);
        catch ME
            if exist('solveTimer', 'var') && ~isempty(solveTimer)
                runResults(k).runtimeFMAM = toc(solveTimer);
            end
            runResults(k).success = false;
            runResults(k).finalLambda = NaN;
            runResults(k).reason = 'exception';
            runResults(k).errorIdentifier = char(string(ME.identifier));
            runResults(k).errorMessage = char(string(ME.message));

            fprintf('Processed scale=%.6g %s: success=0, lambda=NaN, error=%s\n', ...
                scale, fileName, ME.identifier);
        end

        if mod(k, saveEvery) == 0
            numProcessed = k;
            partialRunResults = runResults(1:k);
            numSuccess = nnz([partialRunResults.success]);
            successRate = numSuccess / numProcessed;
            saveFile = fullfile(resultScaleDir, ...
                sprintf('params_inf_sensitivity_summary_%03d.mat', k));
            save(saveFile, 'scale', 'numProcessed', 'numSuccess', ...
                'successRate', 'partialRunResults');
            fprintf('Saved partial summary: %s\n', saveFile);
        end
    end

    numTotal = numel(runResults);
    numSuccess = nnz([runResults.success]);
    successRate = numSuccess / numTotal;

    save(summaryFile, 'scale', 'numTotal', 'numSuccess', 'successRate', 'runResults');

    fprintf('Processed scale=%.6g: %d files, %d succeeded, success rate = %.6f\n', ...
        scale, numTotal, numSuccess, successRate);
    fprintf('Saved summary: %s\n', summaryFile);
end

function dirName = scale_dir_name(scale)
token = strrep(sprintf('%.12g', scale), '.', 'p');
dirName = ['scale_' token];
end
