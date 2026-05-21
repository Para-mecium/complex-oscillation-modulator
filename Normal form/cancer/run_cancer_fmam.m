% Run FMAM AM/FM sweeps for the cancer normal-form example.

scriptDir = fileparts(mfilename('fullpath'));
normalFormDir = fileparts(scriptDir);
repoDir = fileparts(normalFormDir);
addpath(repoDir);
addpath(normalFormDir);
addpath(scriptDir);
addpath(fullfile(repoDir, 'PO_extract'));

baseParameters = [0.35, 5, 1, 2.5, 0.15, 0, 0, 0];
controlledIdx = [6, 7, 8];
initialCondition = [2.1; 2.9];
phaseVariable = struct('name', 'var', 'idx', 1);
truncationOrder = 50;
errBound = 1e-6;
finalErrBound = 1e-12;
amScaleList = 0.8:0.01:1.2;
periodScaleList = 0.7:0.02:1.3;
outputFile = fullfile(scriptDir, 'cancer_fmam_results.mat');

orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'tspan', [0, 1000], ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-8, 'AbsTol', 1e-10), ...
    'minCrossings', 3, ...
    'transientFraction', 0.5);

continuationOptions = struct( ...
    'predictorMode', 'constant', ...
    'initialLambdaStep', 0.01, ...
    'maxLambdaStep', 0.01);

run(fullfile(scriptDir, 'System.m'));
obs = {};
odeFunc = @(t, y, parameter) [ ...
    sys{1}(reshape(y, 1, []), parameter); ...
    sys{2}(reshape(y, 1, []), parameter)];

poResult = extract_periodic_orbit(odeFunc, initialCondition, baseParameters, orbitOptions);
baseState = state(obs, baseParameters, poResult.orbit_t, poResult.orbit_y, truncationOrder, phaseVariable);
baseState.updatePeriod();
baseState.updateVar2();
baseSolverView = fmam_state_ops.solverViewFromState(baseState);
baseDerivedView = fmam_state_ops.derivedViewFromState(baseState);
derivatives = build_symbolic_derivatives(sys, obs, numel(baseParameters));

baseline.params = reshape(baseSolverView.params, 1, []);
baseline.period = baseDerivedView.period;
baseline.varAmp = reshape(baseDerivedView.varAmp, 1, []);

am = struct([]);
solverView = baseSolverView;
for i = 1:numel(amScaleList)
    targetPeriod = baseline.period;
    targetVarAmp = baseline.varAmp .* amScaleList(i);

    itemsPer = struct([]);
    itemsPer(1).prop = 'p_Psi';
    itemsPer(1).idx = 1;
    itemsPer(1).target = targetPeriod / (2 * pi);
    itemsPer(2).prop = 'varAmp';
    itemsPer(2).idx = 1;
    itemsPer(2).target = targetVarAmp(1);
    itemsPer(3).prop = 'varAmp';
    itemsPer(3).idx = 2;
    itemsPer(3).target = targetVarAmp(2);

    task = FMAM_ODE(sys, obs, solverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, ...
        'continuationOptions', continuationOptions, ...
        'verbose', false);
    task.psiUpdateMode = true;
    task.refreshPsiModeReferences();

    task.fit();
    task.step();
    task.errBound = finalErrBound;
    task.fit();

    solverView = task.exportSolverView();
    derivedView = task.exportDerivedView();

    am(i).scale = amScaleList(i);
    am(i).targetPeriod = targetPeriod;
    am(i).targetVarAmp = targetVarAmp;
    am(i).params = reshape(solverView.params, 1, []);
    am(i).finalPeriod = derivedView.period;
    am(i).finalVarAmp = reshape(derivedView.varAmp, 1, []);
    am(i).lambda = task.continuationStatus.lambda;
    am(i).reason = task.continuationStatus.reason;

    fprintf('AM scale %.3g: lambda=%.3g, status=%s\n', ...
        am(i).scale, am(i).lambda, am(i).reason);
end

fm = struct([]);
solverView = baseSolverView;
for i = 1:numel(periodScaleList)
    targetPeriod = baseline.period .* periodScaleList(i);
    targetVarAmp = baseline.varAmp;

    itemsPer = struct([]);
    itemsPer(1).prop = 'p_Psi';
    itemsPer(1).idx = 1;
    itemsPer(1).target = targetPeriod / (2 * pi);
    itemsPer(2).prop = 'varAmp';
    itemsPer(2).idx = 1;
    itemsPer(2).target = targetVarAmp(1);
    itemsPer(3).prop = 'varAmp';
    itemsPer(3).idx = 2;
    itemsPer(3).target = targetVarAmp(2);

    task = FMAM_ODE(sys, obs, baseSolverView, itemsPer, controlledIdx, [], errBound, ...
        'derivatives', derivatives, ...
        'continuationOptions', continuationOptions, ...
        'verbose', false);
    task.psiUpdateMode = true;
    task.refreshPsiModeReferences();

    task.fit();
    task.step();
    task.errBound = finalErrBound;
    task.fit();

    solverView = task.exportSolverView();
    derivedView = task.exportDerivedView();

    fm(i).scale = periodScaleList(i);
    fm(i).targetPeriod = targetPeriod;
    fm(i).targetVarAmp = targetVarAmp;
    fm(i).params = reshape(solverView.params, 1, []);
    fm(i).finalPeriod = derivedView.period;
    fm(i).finalVarAmp = reshape(derivedView.varAmp, 1, []);
    fm(i).lambda = task.continuationStatus.lambda;
    fm(i).reason = task.continuationStatus.reason;

    fprintf('FM scale %.3g: lambda=%.3g, status=%s\n', ...
        fm(i).scale, fm(i).lambda, fm(i).reason);
end

result.baseline = baseline;
result.am = am;
result.fm = fm;

save(outputFile, 'result');
fprintf('Saved cancer FMAM results: %s\n', outputFile);
