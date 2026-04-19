clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
dataDir = fullfile(scriptDir, 'data', 'fig3b');
outputFile = fullfile(dataDir, 'fig3b_base_model_data.mat');

addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));
addpath(genpath(fullfile(repoDir, 'MatCont7p6')));
mkdir(dataDir);

%% Base model for Fig. 3b
initialState = [1; 0];

I1Values = 0.8:0.1:2.0;
ETValues = 0.6:0.2:2.4;

fixedETForI1Slice = 1;
fixedI1ForETSlice = 1;

representativeI1 = 1;
representativeETValues = [0.9, 1.2, 1.5];

%% Periodic-orbit extraction settings
singleTimeSpan = 1600;
maxWindows = 3;
solverTol = struct('RelTol', 1e-6, 'AbsTol', 1e-9);

orbitOptions = struct();
orbitOptions.solver_name = 'ode45';
orbitOptions.single_timespan = singleTimeSpan;
orbitOptions.max_windows = maxWindows;
orbitOptions.event = 1;
orbitOptions.solver_tol = solverTol;
orbitOptions.minCrossings = 6;
orbitOptions.transientFraction = 0.5;
orbitOptions.samplesPerCycle = 400;
orbitOptions.extractNumPoints = 500;

%% Scan the (I1, ET) grid
periodGrid = NaN(numel(ETValues), numel(I1Values));
amplitudeGrid = NaN(numel(ETValues), numel(I1Values));
gridHasOrbit = false(numel(ETValues), numel(I1Values));

for i = 1:numel(I1Values)
    for j = 1:numel(ETValues)
        params = [I1Values(i), ETValues(j)];
        fprintf('Grid point %d/%d, %d/%d: I1 = %.2f, ET = %.2f\n', ...
            i, numel(I1Values), j, numel(ETValues), params(1), params(2));

        result = flexmod_forward_orbit(params, initialState, struct( ...
            'systemName', 'base', ...
            'poOptions', orbitOptions));

        gridHasOrbit(j, i) = result.success;
        if result.success
            periodGrid(j, i) = result.features.period;
            amplitudeGrid(j, i) = result.features.state.amplitude(2);
        end
    end
end

%% Sweep I1 with ET fixed at 1
I1SliceValues = I1Values(:);
I1SlicePeriods = NaN(numel(I1SliceValues), 1);
I1SliceAmplitudes = NaN(numel(I1SliceValues), 1);
I1SliceHasOrbit = false(numel(I1SliceValues), 1);

for i = 1:numel(I1SliceValues)
    params = [I1SliceValues(i), fixedETForI1Slice];
    result = flexmod_forward_orbit(params, initialState, struct( ...
        'systemName', 'base', ...
        'poOptions', orbitOptions));

    I1SliceHasOrbit(i) = result.success;
    if result.success
        I1SlicePeriods(i) = result.features.period;
        I1SliceAmplitudes(i) = result.features.state.amplitude(2);
    end
end

%% Sweep ET with I1 fixed at 1
ETSliceValues = ETValues(:);
ETSlicePeriods = NaN(numel(ETSliceValues), 1);
ETSliceAmplitudes = NaN(numel(ETSliceValues), 1);
ETSliceHasOrbit = false(numel(ETSliceValues), 1);

for i = 1:numel(ETSliceValues)
    params = [fixedI1ForETSlice, ETSliceValues(i)];
    result = flexmod_forward_orbit(params, initialState, struct( ...
        'systemName', 'base', ...
        'poOptions', orbitOptions));

    ETSliceHasOrbit(i) = result.success;
    if result.success
        ETSlicePeriods(i) = result.features.period;
        ETSliceAmplitudes(i) = result.features.state.amplitude(2);
    end
end

%% Representative time traces at I1 = 1
representativeOrbits = cell(1, numel(representativeETValues));

for i = 1:numel(representativeETValues)
    params = [representativeI1, representativeETValues(i)];
    result = flexmod_forward_orbit(params, initialState, struct( ...
        'systemName', 'base', ...
        'poOptions', orbitOptions, ...
        'shiftToYMax', true));

    representativeOrbits{i} = result.orbit;
end

%% Save data
save(outputFile, ...
    'I1Values', 'ETValues', ...
    'periodGrid', 'amplitudeGrid', 'gridHasOrbit', ...
    'I1SliceValues', 'fixedETForI1Slice', 'I1SlicePeriods', 'I1SliceAmplitudes', 'I1SliceHasOrbit', ...
    'ETSliceValues', 'fixedI1ForETSlice', 'ETSlicePeriods', 'ETSliceAmplitudes', 'ETSliceHasOrbit', ...
    'representativeI1', 'representativeETValues', 'representativeOrbits');

fprintf('Saved data: %s\n', outputFile);
