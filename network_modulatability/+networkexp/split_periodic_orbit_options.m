function [solverOptions, detectionOptions] = split_periodic_orbit_options(options)
%SPLIT_PERIODIC_ORBIT_OPTIONS Separate ODE solver options from orbit-detection options.

if nargin < 1 || isempty(options)
    solverOptions = struct();
    detectionOptions = struct();
    return
end
if ~isstruct(options)
    error('networkexp:InvalidPeriodicOrbitOptions', ...
        'Periodic-orbit options must be empty or a struct.');
end

detectionFieldNames = { ...
    'observableFcn', 'autoSection', 'sectionFcn', 'sectionLevel', ...
    'sectionDirection', 'transientFraction', 'transientTime', ...
    'minCrossings', 'consecutiveCycles', 'poincareTol', 'periodTol', ...
    'amplitudeTol', 'poincareTolLoose', 'periodTolLoose', ...
    'amplitudeTolLoose', 'nonoscAmpTol', 'nonoscStdTol', ...
    'minPeakProminence', 'minPeakDistance', 'samplesPerCycle', ...
    'extractNumPoints', 'computeQualityDiagnostics', 'minEventSeparation', 'event'};

solverOptions = options;
detectionOptions = struct();
for idxField = 1:numel(detectionFieldNames)
    fieldName = detectionFieldNames{idxField};
    if isfield(solverOptions, fieldName)
        detectionOptions.(fieldName) = solverOptions.(fieldName); %#ok<AGROW>
        solverOptions = rmfield(solverOptions, fieldName);
    end
end
end
