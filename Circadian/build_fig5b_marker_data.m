clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'data', 'common', 'initData.mat');
markerDataDir = fullfile(scriptDir, 'data', 'fig5b', 'markers');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(markerDataDir);

%% Fig. 5b marker settings
targetPeriod = 24.0;
targetAmplitudes = [0.01, 0.02, 0.04];
markSeeds = [5e-5, 0.06; ...
             5e-5, 0.06; ...
             5e-5, 0.06];

obs = {@(variable) variable(:, 2) + variable(:, 3)};
M = 50;
PV.name = 'obs';
PV.idx = 1;

controlledIdx = [1, 2];
errBound = 1e-6;

continuationOptions = struct( ...
    'initialLambdaStep', 0.05, ...
    'predictorMode', 'constant', ...
    'conditionStopEnabled', true, ...
    'conditionStopRcond', 1e-9);

%% Build system and derivatives
initData = load(initDataFile);
baseY0 = initData.TS{2}(1, :).';

run(fullfile(scriptDir, 'System.m'));
derivatives = build_symbolic_derivatives(sys, obs, 2);

%% Compute marker points and save data
for i = 1:numel(targetAmplitudes)
    targetAmplitude = targetAmplitudes(i);
    startParams = markSeeds(i, :);
    saveFile = fullfile(markerDataDir, ...
        sprintf('fig5b_marker_T%s_A%s.mat', period_tag(targetPeriod), amplitude_tag(targetAmplitude)));

    startOrbit = circadian_forward_orbit(startParams, baseY0, struct());
    if ~startOrbit.success
        error('circadian_refactor:Fig5bMarkerSeedOrbitGenerationFailed', ...
            'Periodic-orbit extraction failed at marker seed %d.', i);
    end

    seedState = state(obs, startParams, startOrbit.orbit.t, startOrbit.orbit.y, M, PV);
    seedState.updatePeriod();
    seedState.updateVar2();
    seedSolverView = fmam_state_ops.solverViewFromState(seedState);

    itemsPer = struct([]);
    itemsPer(1).prop = 'p_Psi';
    itemsPer(1).idx = 1;
    itemsPer(1).target = targetPeriod / (2 * pi);

    itemsPer(2).prop = 'obsAmp';
    itemsPer(2).idx = 1;
    itemsPer(2).target = targetAmplitude;

    markerTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    markerTask.psiUpdateMode = true;
    markerTask.refreshPsiModeReferences();
    markerTask.needLog = false;

    markerTask.fit()
    markerTask.step()
    markerTask.errBound = 1e-12;
    markerTask.fit()

    markerDerivedView = markerTask.exportDerivedView();
    Parameters = reshape(markerTask.exportSolverView().params, 1, []);

    markerOrbit = circadian_forward_orbit(Parameters, markerDerivedView.TS_var(1, :).', struct());
    if ~markerOrbit.success
        error('circadian_refactor:Fig5bMarkerOrbitGenerationFailed', ...
            'Periodic-orbit extraction failed at marker point %d.', i);
    end

    TS = {markerOrbit.orbit.t, markerOrbit.orbit.y};
    period = markerOrbit.features.period;
    obsAmp = reshape(markerOrbit.features.observable.amplitude, 1, []);
    obsMax = reshape(markerOrbit.features.observable.max, 1, []);
    obsMin = reshape(markerOrbit.features.observable.min, 1, []);

    save(saveFile, ...
        'Parameters', ...
        'TS', ...
        'period', ...
        'obsAmp', ...
        'obsMax', ...
        'obsMin');

    fprintf('Saved marker data: %s\n', saveFile);
end

%%
function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end

function tag = amplitude_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end
