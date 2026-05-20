clear
clc
%%
selectedNoiseLevels = [0.01 0.02 0.05];
selectedFiles = 1:100;
%%
scriptDir = fileparts(mfilename('fullpath'));
parameterInferenceDir = fileparts(scriptDir);
circuitDir = fileparts(parameterInferenceDir);
repoDir = fileparts(circuitDir);
inputDir = fullfile(scriptDir, 'noisy_init_data_files');
outputDir = fullfile(scriptDir, 'modulation_results');

addpath(repoDir, '-begin');
addpath(circuitDir, '-begin');

initData = load(fullfile(parameterInferenceDir, 'initData_ODE.mat'));
run(fullfile(circuitDir, 'System.m'));

obs = [];
params = initData.Parameters;
t = initData.TS{1};
TS_var = initData.TS{2};
M = 50;
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

PV.name = 'var';
PV.idx = 1;
continuationOptions = struct('initialLambdaStep', 0.01, 'predictorMode', 'constant');
items_controlled = [1 4 5 6];
errBound = 1e-6;

poExtractDir = fullfile(repoDir, 'PO_extract');
addpath(poExtractDir, '-begin');

odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter)];

orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

for noiseIdx = 1:numel(selectedNoiseLevels)
    noiseLevel = selectedNoiseLevels(noiseIdx);
    noiseLevelDirName = ['noise_level_' strrep(num2str(noiseLevel), '.', 'p')];
    sourceNoiseDir = fullfile(inputDir, noiseLevelDirName);
    targetNoiseDir = fullfile(outputDir, noiseLevelDirName);
    summaryFile = fullfile(targetNoiseDir, 'modulation_summary.mat');
    sampleSuccess = false(1, numel(selectedFiles));
    identifiedParameters = cell(1, numel(selectedFiles));

    if ~isfolder(targetNoiseDir)
        mkdir(targetNoiseDir);
    end

    for fileIdx = 1:numel(selectedFiles)
        seedIdx = selectedFiles(fileIdx);
        fileName = sprintf('initData_%03d.mat', seedIdx);
        targetData = load(fullfile(sourceNoiseDir, fileName));

        State = state(obs, params, t, TS_var, M, PV);
        StateView = fmam_state_ops.solverViewFromState(State);

        items_per = struct;
        items_per(1).prop = 'p_Psi';
        items_per(1).idx = 1;
        items_per(1).target = targetData.period / (2 * pi);
        items_per(2).prop = 'varAmp';
        items_per(2).idx = 1;
        items_per(2).target = targetData.varAmp(1);
        items_per(3).prop = 'varAmp';
        items_per(3).idx = 2;
        items_per(3).target = targetData.varAmp(2);
        items_per(4).prop = 'varAmp';
        items_per(4).idx = 3;
        items_per(4).target = targetData.varAmp(3);

        Modtask = FMAM_ODE(sys, obs, StateView, items_per, items_controlled, [], errBound, ...
            'derivatives', derivatives, 'continuationOptions', continuationOptions);
        Modtask.psiUpdateMode = true;
        Modtask.refreshPsiModeReferences();
        Modtask.needLog = false;
        Modtask.verbose = false;

        Modtask.fit();
        Modtask.step();

        finalLambda = Modtask.continuationStatus.lambda;
        finalView = Modtask.exportSolverView();

        success = abs(finalLambda - 1) <= 1e-12;
        Parameters = reshape(finalView.params, 1, []);
        TS = [];
        reason = [];

        if success
            StateDerived = Modtask.exportDerivedView();
            y0 = StateDerived.TS_var(1, :).';
            searchWindow = max(10 * StateDerived.period, StateDerived.t(end));
            orbitOptions.tspan = [0, searchWindow];

            poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
            if poResult.has_orbit
                TS = {poResult.orbit_t, poResult.orbit_y};
            else
                success = false;
                reason = 'orbit_extraction_failed';
            end
        else
            reason = 'continuation_failed';
        end

        saveName = ['modulated_to_' fileName];
        saveFile = fullfile(targetNoiseDir, saveName);
        noiseLevel = targetData.noiseLevel;
        seed = targetData.seed;
        save(saveFile, 'success', 'Parameters', 'TS', 'noiseLevel', 'seed', 'reason');

        sampleSuccess(fileIdx) = success;
        identifiedParameters{fileIdx} = Parameters;

        msg = sprintf('Processed %s: success=%d, lambda=%.6g, output=%s', ...
            fileName, success, finalLambda, saveFile);
        fprintf('%s\n', msg);
    end

    successRate = nnz(sampleSuccess) / numel(sampleSuccess);
    save(summaryFile, 'successRate', 'sampleSuccess', 'identifiedParameters', ...
        'noiseLevel', 'selectedFiles');
end
