clear
clc
%% Extended Data Fig. 3b initial periodic orbit
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
saveFile = fullfile(scriptDir, 'initData.mat');

Parameters = [30.5, 183, 0.1, 0.1, 3.7, 3.7, 0.3, 0.2, 3.8, 326, 185, 3, 4.8];
y0 = [23.3896986043691; 222.452142059269; 207.783959243676; 216.598138352976];
orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 400], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

addpath(fullfile(repoDir, 'PO_extract'));

sys = build_longevity_system();
odeFunc = @(t, y, parameter) ode_rhs_from_sys(sys, y, parameter);
poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
if ~poResult.has_orbit
    error('build_ed_fig3b_init_data:PeriodicOrbitGenerationFailed', ...
        'Periodic-orbit extraction did not return an orbit (%s).', ...
        poResult.message);
end

orbitForFeatures = struct( ...
    't', poResult.orbit_t(:), ...
    'y', poResult.orbit_y, ...
    'period', poResult.period);
poFeatures = evaluate_orbit_features(orbitForFeatures, [], [], struct());

TS = {poResult.orbit_t, poResult.orbit_y};
period = poFeatures.period;
varAmp = reshape(poFeatures.state.amplitude, 1, []);
varMax = reshape(poFeatures.state.max, 1, []);
varMin = reshape(poFeatures.state.min, 1, []);

save(saveFile, 'TS', 'Parameters', 'period', 'varAmp', 'varMax', 'varMin');
fprintf('Saved init data: %s\n', saveFile);

function sys = build_longevity_system()
sys = cell(1, 4);
sys{1} = @dmSdt;
sys{2} = @dmHdt;
sys{3} = @dSdt;
sys{4} = @dHdt;
end

function dydt = ode_rhs_from_sys(sys, y, parameter)
yRow = reshape(y, 1, []);
dydt = zeros(numel(sys), 1);
for i = 1:numel(sys)
    dydt(i) = sys{i}(yRow, parameter);
end
end

function output = dmSdt(TS_variable, parameter)
alphaS = parameter(1);
alphaS0 = parameter(3);
deltam = parameter(7);
KH = parameter(10);
n1 = parameter(12);

mS = TS_variable(:, 1);
H = TS_variable(:, 4);

output = alphaS0 + alphaS * H.^n1 ./ (KH^n1 + H.^n1) - deltam * mS;
end

function output = dmHdt(TS_variable, parameter)
alphaH = parameter(2);
alphaH0 = parameter(4);
deltam = parameter(7);
KS = parameter(11);
n2 = parameter(13);

mH = TS_variable(:, 2);
S = TS_variable(:, 3);

output = alphaH0 + alphaH * KS^n2 ./ (KS^n2 + S.^n2) - deltam * mH;
end

function output = dSdt(TS_variable, parameter)
betaS = parameter(5);
deltaS = parameter(8);

mS = TS_variable(:, 1);
S = TS_variable(:, 3);

output = betaS * mS - deltaS * S;
end

function output = dHdt(TS_variable, parameter)
betaH = parameter(6);
deltaH = parameter(9);

mH = TS_variable(:, 2);
H = TS_variable(:, 4);

output = betaH * mH - deltaH * H;
end
