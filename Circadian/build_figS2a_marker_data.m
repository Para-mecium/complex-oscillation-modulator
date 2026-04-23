clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
initDataFile = fullfile(scriptDir, 'data', 'common', 'initData.mat');
markerDataDir = fullfile(scriptDir, 'data', 'figS2a', 'markers');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
mkdir(markerDataDir);

%% Fig. S2a marker settings
targetAmplitude = 0.025;
targetPeriods = [23.5, 24.0, 24.5];
markSeeds = [1.23601799538e-4, 0.0640753316922; ...
             7.86535122639e-5, 0.0487912547039; ...
             5.10865154102e-5, 0.0404189029782];

%% FMAM and observable settings
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
for i = 1:numel(targetPeriods)
    targetPeriod = targetPeriods(i);
    startParams = markSeeds(i, :);
    saveFile = fullfile(markerDataDir, ...
        sprintf('figS2a_marker_A%s_T%s.mat', amplitude_tag(targetAmplitude), period_tag(targetPeriod)));

    startOrbit = circadian_forward_orbit(startParams, baseY0, struct());
    if ~startOrbit.success
        error('circadian_refactor:FigS2aMarkerSeedOrbitGenerationFailed', ...
            'Periodic-orbit extraction failed at marker seed %d.', i);
    end

    seedState = state(obs, startParams, startOrbit.orbit.t, startOrbit.orbit.y, M, PV);
    seedState.updatePeriod();
    seedState.updateVar2();
    seedSolverView = fmam_state_ops.solverViewFromState(seedState);

    itemsPer = struct([]);
    itemsPer(1).prop = 'obsAmp';
    itemsPer(1).idx = 1;
    itemsPer(1).target = targetAmplitude;

    itemsPer(2).prop = 'p_Psi';
    itemsPer(2).idx = 1;
    itemsPer(2).target = targetPeriod / (2 * pi);

    markerTask = FMAM_ODE(sys, obs, seedSolverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, 'continuationOptions', continuationOptions);
    markerTask.psiUpdateMode = true;
    markerTask.refreshPsiModeReferences();
    markerTask.needLog = false;

    markerTask.fit()
    markerTask.step()
    markerTask.errBound = 1e-12;
    markerTask.fit()

    Parameters = reshape(markerTask.exportSolverView().params, 1, []);
    markerDerivedView = markerTask.exportDerivedView();

    markerOrbit = circadian_forward_orbit(Parameters, markerDerivedView.TS_var(1, :).', struct());
    if ~markerOrbit.success
        error('circadian_refactor:FigS2aMarkerOrbitGenerationFailed', ...
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
function tag = amplitude_tag(value)
tag = strrep(sprintf('%.3f', value), '.', 'p');
end

function tag = period_tag(value)
tag = strrep(sprintf('%.1f', value), '.', 'p');
end
