function options = build_fhn_orbit_detection_options(baseOptions, searchWindow, weights)
%BUILD_FHN_ORBIT_DETECTION_OPTIONS Configure network-level periodic-orbit detection.

if nargin < 1 || isempty(baseOptions)
    baseOptions = struct();
end
if nargin < 2
    searchWindow = [];
end
if nargin < 3
    weights = [];
end

options = baseOptions;

if isempty(weights)
    options.observableFcn = @mean_voltage_observable;
else
    weights = weights(:);
    weights = weights / sum(weights);
    options.observableFcn = @(t, y, params) weighted_voltage_observable(t, y, params, weights); %#ok<INUSD>
end

options.autoSection = false;
options.sectionLevel = 0.3;
options.sectionDirection = +1;
options.poincareTol = 0.15;
if ~isempty(searchWindow)
    options.transientTime = 0.5 * searchWindow;
end
end

function value = mean_voltage_observable(~, y, params)
numNodes = size(params{6}, 1);
value = mean(y(1:numNodes));
end

function value = weighted_voltage_observable(~, y, ~, weights)
value = weights.' * y(1:numel(weights));
end
