clear
clc

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(scriptDir, 'data', 'bifurcation');
outputFile = fullfile(dataDir, 'circadian_bifurcation_line.mat');

mkdir(dataDir);

%% Bifurcation-line settings
beta = 0.1572;
nHopfPoints = 3000;
xRange = [0, 1.5e-4];
yRange = [0, 0.09];
equilibriumX = linspace(0, 0.12, nHopfPoints).';

%% Compute Hopf curve
AT = equilibriumX .* (equilibriumX + 7) ./ (8 * (1 - equilibriumX).^2);
Kd = AT .* (1 - 2 * equilibriumX) - 7 * equilibriumX / 8;
feedbackSlope = -8 * ones(size(equilibriumX));
hopfFrequency = sqrt(3) * beta * ones(size(equilibriumX));

positiveMask = Kd >= 0 & AT >= 0;
hopfCurve = struct();
hopfCurve.Kd = Kd(positiveMask);
hopfCurve.AT = AT(positiveMask);
hopfCurve.equilibrium = equilibriumX(positiveMask);
hopfCurve.feedbackSlope = feedbackSlope(positiveMask);
hopfCurve.frequency = hopfFrequency(positiveMask);

visibleMask = positiveMask & ...
    Kd >= xRange(1) & Kd <= xRange(2) & ...
    AT >= yRange(1) & AT <= yRange(2);

visibleHopfCurve = struct();
visibleHopfCurve.Kd = Kd(visibleMask);
visibleHopfCurve.AT = AT(visibleMask);
visibleHopfCurve.equilibrium = equilibriumX(visibleMask);
visibleHopfCurve.feedbackSlope = feedbackSlope(visibleMask);
visibleHopfCurve.frequency = hopfFrequency(visibleMask);

save(outputFile, ...
    'beta', ...
    'xRange', ...
    'yRange', ...
    'hopfCurve', ...
    'visibleHopfCurve');

fprintf('Saved bifurcation data: %s\n', outputFile);
fprintf('Visible Hopf points: %d\n', numel(visibleHopfCurve.Kd));
