function tests = test_utils_state_ops
%TEST_UTILS_STATE_OPS Regression coverage for canonical FMAM state helpers.

    tests = functiontests(localfunctions);
end

function testSeriesReconstructionMatchesStateProperties(testCase)
    stat = make_reference_state();
    phaseSampleCount = stat.discretization.reconstruction.phaseSampleCount;

    TS_var = fmam_state_ops.evaluateTrigSeries(stat.phi, stat.p_var, stat.q_var);
    Psi = fmam_state_ops.evaluateTrigSeries(stat.phi, stat.p_Psi, stat.q_Psi);
    TS_obs = fmam_state_ops.getObs(stat.obs, TS_var);
    dt = fmam_state_ops.timeIncrementsFromCoefficients(stat.p_Psi, stat.q_Psi, phaseSampleCount);
    t = [0; cumsum(dt(1:end-1))];

    verifyEqual(testCase, TS_var, stat.TS_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, Psi, stat.Psi, 'AbsTol', 1e-10);
    verifyEqual(testCase, TS_obs, stat.TS_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, t, stat.t, 'AbsTol', 1e-10);
end

function testDerivedStateBuildersMatchStateStorage(testCase)
    stat = make_reference_state();
    extremaSearch = stat.extremaSearch;
    discretization = stat.discretization;
    phaseSampleCount = discretization.reconstruction.phaseSampleCount;

    varDerived = fmam_state_ops.buildVariableDerivedState( ...
        stat.p_var, stat.q_var, stat.p_Psi, stat.q_Psi, ...
        phaseSampleCount, extremaSearch);
    obsDerived = fmam_state_ops.buildObservableDerivedState( ...
        stat.obs, stat.p_var, stat.q_var, stat.p_Psi, stat.q_Psi, ...
        discretization, extremaSearch);
    [pOrigin, qOrigin] = fmam_state_ops.reprojectEqualTimeFourier( ...
        stat.t / (stat.period / (2 * pi)), stat.TS_var, stat.truncationOrder);

    verifyEqual(testCase, varDerived.varPhiMax, stat.varPhiMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, varDerived.varPhiMin, stat.varPhiMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, varDerived.varAmp, stat.varAmp, 'AbsTol', 1e-10);
    verifyEqual(testCase, varDerived.varMax, stat.varMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, varDerived.varMin, stat.varMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, varDerived.varPhase, stat.varPhase, 'AbsTol', 1e-10);

    verifyEqual(testCase, obsDerived.p_obs, stat.p_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.q_obs, stat.q_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.obsPhiMax, stat.obsPhiMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.obsPhiMin, stat.obsPhiMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.obsAmp, stat.obsAmp, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.obsMax, stat.obsMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.obsMin, stat.obsMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, obsDerived.obsPhase, stat.obsPhase, 'AbsTol', 1e-10);

    verifyEqual(testCase, pOrigin, stat.p_var_origin, 'AbsTol', 1e-10);
    verifyEqual(testCase, qOrigin, stat.q_var_origin, 'AbsTol', 1e-10);
end

function testSolverViewBuilderMatchesStateSolverStorage(testCase)
    [obs,t,x] = reference_trajectory();
    stat = state(obs, [1 2], t, x, 5, struct('name', 'var', 'idx', 1));

    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        obs,[1 2],t,x,5,struct('name', 'var', 'idx', 1), ...
        stat.discretization, stat.extremaSearch);

    verifyEqual(testCase, solverView.params, stat.params, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.p_Psi, stat.p_Psi, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.q_Psi, stat.q_Psi, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.p_var, stat.p_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.q_var, stat.q_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.varPhiMax, stat.varPhiMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.varPhiMin, stat.varPhiMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.obsPhiMax, stat.obsPhiMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.obsPhiMin, stat.obsPhiMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, size(solverView.p_Psi), [stat.truncationOrder, 1]);
    verifyEqual(testCase, size(solverView.q_Psi), [stat.truncationOrder - 1, 1]);
end

function testDerivedViewBuilderMatchesStateStorage(testCase)
    stat = make_reference_state();
    solverView = struct( ...
        'params', stat.params, ...
        'p_Psi', stat.p_Psi, ...
        'q_Psi', stat.q_Psi, ...
        'p_var', stat.p_var, ...
        'q_var', stat.q_var, ...
        'varPhiMax', stat.varPhiMax, ...
        'varPhiMin', stat.varPhiMin, ...
        'obsPhiMax', stat.obsPhiMax, ...
        'obsPhiMin', stat.obsPhiMin, ...
        'PV', stat.PV);

    derived = fmam_state_ops.buildDerivedView( ...
        stat.obs,solverView,stat.discretization);

    verifyEqual(testCase, derived.Psi, stat.Psi, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.t, stat.t, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.TS_var, stat.TS_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.TS_obs, stat.TS_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.period, stat.period, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.p_var_origin, stat.p_var_origin, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.q_var_origin, stat.q_var_origin, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.p_obs, stat.p_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.q_obs, stat.q_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.varAmp, stat.varAmp, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.varMax, stat.varMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.varMin, stat.varMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.varPhase, stat.varPhase, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.obsAmp, stat.obsAmp, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.obsMax, stat.obsMax, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.obsMin, stat.obsMin, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.obsPhase, stat.obsPhase, 'AbsTol', 1e-10);
end

function testStateViewExportHelpersMatchStateStorage(testCase)
    stat = make_reference_state();

    solverView = fmam_state_ops.solverViewFromState(stat);
    derived = fmam_state_ops.derivedViewFromState(stat);
    snapshot = fmam_state_ops.buildStateSnapshotFromViews(solverView, derived);

    verifyEqual(testCase, solverView.params, stat.params, 'AbsTol', 1e-10);
    verifyEqual(testCase, solverView.p_var, stat.p_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.TS_var, stat.TS_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.period, stat.period, 'AbsTol', 1e-10);
    verifyEqual(testCase, snapshot.p_var_origin, stat.p_var_origin, 'AbsTol', 1e-10);
    verifyEqual(testCase, snapshot.obsPhase, stat.obsPhase, 'AbsTol', 1e-10);
end

function testDerivedViewFromSnapshotReconstructsCanonicalSeries(testCase)
    stat = make_reference_state();
    phaseSampleCount = stat.discretization.reconstruction.phaseSampleCount;
    snapshot = fmam_state_ops.buildStateSnapshotFromViews( ...
        fmam_state_ops.solverViewFromState(stat), ...
        fmam_state_ops.derivedViewFromState(stat));

    derived = fmam_state_ops.derivedViewFromSnapshot(snapshot, phaseSampleCount);

    verifyEqual(testCase, derived.Psi, stat.Psi, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.t, stat.t, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.TS_var, stat.TS_var, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.TS_obs, stat.TS_obs, 'AbsTol', 1e-10);
    verifyEqual(testCase, derived.period, stat.period, 'AbsTol', 1e-10);
end

function testValidateTrajectoryInputsAcceptsCanonicalTrajectory(testCase)
    [obs,t,x] = reference_trajectory();

    verifyWarningFree(testCase, @() fmam_state_ops.validateTrajectoryInputs( ...
        obs, [1 2], t, x, 5, struct('name', 'var', 'idx', 1)));
end

function testHelpersAcceptUnifiedDiscretizationStruct(testCase)
    stat = make_reference_state();
    discretization = fmam_state_defaults.defaultDiscretization();
    extremaSearch = fmam_state_ops.defaultExtremaSearchSettings();
    discretization.reconstruction.timeResampleCount = 4096;
    discretization.reconstruction.phaseSampleCount = 64;

    solverView = fmam_state_ops.buildSolverViewFromTrajectory( ...
        stat.obs, stat.params, stat.t, stat.TS_var, stat.truncationOrder, stat.PV, ...
        discretization, extremaSearch);
    derived = fmam_state_ops.buildDerivedView(stat.obs, solverView, discretization);

    verifyEqual(testCase, size(derived.TS_var, 1), discretization.reconstruction.phaseSampleCount);
    verifyEqual(testCase, size(derived.TS_obs, 1), discretization.reconstruction.phaseSampleCount);
    verifyEqual(testCase, solverView.PV, stat.PV);
end

function stat = make_reference_state()
    [obs,t,x] = reference_trajectory();
    stat = state(obs, [1 2], t, x, 5, struct('name', 'var', 'idx', 1));
end

function [obs,t,x] = reference_trajectory()
    obs = {@observable_mix};
    t = linspace(0, 2 * pi, 1001).';
    t(end) = [];
    x = [cos(t), 0.5 * sin(t) + 0.2 * cos(2 * t)];
end

function z = observable_mix(y)
    z = y(:,1) - 0.5 * y(:,2);
end
