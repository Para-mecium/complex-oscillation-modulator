clear
clc

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
addpath(repoDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));

Parameters0 = [2, 300, 0, 0.2, 1, 1];
y0 = [1; 0; 0; 0; 0; 0];
M = 75;
PV.name = 'var';
PV.idx = 4;
periodTarget = 50;
phaseTargets = [30, 35];
errBound = 1e-6;

run(fullfile(scriptDir, 'System.m'));
obs = {};
derivatives = build_symbolic_derivatives(sys, obs, numel(Parameters0));

baseOrbit = repressilator_forward_orbit(Parameters0, y0, struct());
TS = {baseOrbit.orbit.t, baseOrbit.orbit.y};
Parameters = Parameters0;
period = baseOrbit.features.period;
save(fullfile(scriptDir, 'initData_repressilator.mat'), 'Parameters', 'TS', 'period');

baseState = state(obs, Parameters0, baseOrbit.orbit.t, baseOrbit.orbit.y, M, PV);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);

items = struct([]);
items(1).prop = 'p_Psi';
items(1).idx = 1;
items(1).target = periodTarget / (2 * pi);

periodTask = FMAM_ODE(sys, obs, baseSolverView, items, 4, [], errBound, ...
    'derivatives', derivatives, 'verbose', false);
periodTask.psiUpdateMode = true;
periodTask.refreshPsiModeReferences();
periodTask.fit();
periodTask.step();
periodTask.errBound = 1e-12;
periodTask.fit();

periodSolverView = periodTask.exportSolverView();
periodDerivedView = periodTask.exportDerivedView();
periodParameters = reshape(periodSolverView.params, 1, []);
periodOrbit = repressilator_forward_orbit(periodParameters, periodDerivedView.TS_var(1, :).', struct());

periodData = struct();
periodData.Parameters = periodParameters;
periodData.TS = {periodOrbit.orbit.t, periodOrbit.orbit.y};
periodData.period = periodOrbit.features.period;
periodData.phase12 = phase_difference(periodOrbit.orbit.t, periodOrbit.orbit.y(:, 4), periodOrbit.orbit.y(:, 5));
save(fullfile(scriptDir, 'period_target_50.mat'), '-struct', 'periodData');

phaseResults = struct([]);
for i = 1:numel(phaseTargets)
    items = struct([]);
    items(1).prop = 'p_Psi';
    items(1).idx = 1;
    items(1).target = periodTarget / (2 * pi);
    items(2).prop = 'varPhase';
    items(2).idx = [4, 5];
    items(2).target = phaseTargets(i);

    phaseTask = FMAM_ODE(sys, obs, periodSolverView, items, [4, 5], [0.1, 0.1], errBound, ...
        'derivatives', derivatives, 'verbose', false);
    phaseTask.psiUpdateMode = true;
    phaseTask.refreshPsiModeReferences();
    phaseTask.fit();
    phaseTask.step();
    phaseTask.errBound = 1e-12;
    phaseTask.fit();

    solverView = phaseTask.exportSolverView();
    derivedView = phaseTask.exportDerivedView();
    finalParameters = reshape(solverView.params, 1, []);
    finalOrbit = repressilator_forward_orbit(finalParameters, derivedView.TS_var(1, :).', struct());

    phaseResults(i).targetPhase = phaseTargets(i);
    phaseResults(i).Parameters = finalParameters;
    phaseResults(i).TS = {finalOrbit.orbit.t, finalOrbit.orbit.y};
    phaseResults(i).period = finalOrbit.features.period;
    phaseResults(i).phase12 = phase_difference(finalOrbit.orbit.t, finalOrbit.orbit.y(:, 4), finalOrbit.orbit.y(:, 5));
    phaseResults(i).lambda = phaseTask.continuationStatus.lambda;
    phaseResults(i).reason = phaseTask.continuationStatus.reason;

    Parameters = finalParameters;
    TS = phaseResults(i).TS;
    period = phaseResults(i).period;
    phase12 = phaseResults(i).phase12;
    targetPhase = phaseTargets(i);
    save(fullfile(scriptDir, sprintf('phase_target_%d.mat', phaseTargets(i))), ...
        'Parameters', 'TS', 'period', 'phase12', 'targetPhase');
end

save(fullfile(scriptDir, 'figS17_repressilator_data.mat'), ...
    'Parameters0', 'periodTarget', 'phaseTargets', 'periodData', 'phaseResults');

function phase12 = phase_difference(t, p1, p2)
    tq = linspace(t(1), t(end), 20000).';
    p1q = interp1(t, p1, tq, 'spline');
    p2q = interp1(t, p2, tq, 'spline');
    [~, idx1] = max(p1q);
    [~, idx2] = max(p2q);
    phase12 = tq(idx2) - tq(idx1);
    if phase12 < 0
        phase12 = phase12 + t(end);
    end
end
