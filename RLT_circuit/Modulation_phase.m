%% FMAM
% This block of codes performs a tri-modulation task.
% We modulate the period, the amplitude of V1 and the amplitude of V2
% through varying Rc, beta1 and beta2
clear

% load data
System
obs = [];

% load("T_A1_A2_A3.mat")
% load("initData (ODE, Imax = 2.97mA).mat")
load("learnedData_ODE.mat")

params = Parameters;

t = TS{1};
TS_var = TS{2};
TS_obs = [];
M = 75;
cd ../
derivatives = build_symbolic_derivatives(sys, obs, numel(params));

% Set Primary variable
PV.name = 'var';
PV.idx = 1;
State = state(obs,params,t,TS_var,M,PV);
StateView = fmam_state_ops.solverViewFromState(State);
%%
% Set modualtion target
period_target = 30;
amplitude1_target = State.params(3)*10;
items_per = struct;

items_per(1).prop = 'p_Psi';
items_per(1).idx = 1;
items_per(1).target = period_target/(2*pi);

items_per(2).prop = 'p_var';
items_per(2).idx = [2,1];
items_per(2).target = amplitude1_target;

% Set control parameters
items_controlled = [1 4];
% Set stepsize
max_stepsize = abs([(period_target-State.period)/(2*pi) amplitude1_target-State.varAmp(1)])/100;
errBound = 1e-6;

Modtask = FMAM_ODE(sys,[],StateView,items_per,items_controlled,[],errBound, ...
    'derivatives', derivatives);
Modtask.psiUpdateMode = true;
Modtask.refreshPsiModeReferences();
Modtask.needLog = true;
%%
tic
Modtask.fit()
Modtask.step()
elapsedTime = toc;  % 
disp(['Computing time: ', num2str(elapsedTime), ' seconds']);

StateView = Modtask.exportSolverView();
StateDerived = Modtask.exportDerivedView();

%%
addpath(fullfile(pwd, 'PO_extract'), '-begin');
odeFunc = @(t, y, parameter) [ ...
    sys{1}(y(:).', parameter); ...
    sys{2}(y(:).', parameter); ...
    sys{3}(y(:).', parameter)];
orbitOptions = struct( ...
    'solver_name', 'ode45', ...
    'event', 1, ...
    'solver_tol', struct('RelTol', 1e-6, 'AbsTol', 1e-9), ...
    'minCrossings', 3, ...
    'transientFraction', 0);

baseStateView = StateView;
phase12_targets = [20 25];

figure
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for targetIdx = 1:numel(phase12_targets)
    % Set modualtion target
    phase12_target = phase12_targets(targetIdx);
    items_per = struct;

    items_per(1).prop = 'p_Psi';
    items_per(1).idx = 1;
    items_per(1).target = period_target/(2*pi);

    items_per(2).prop = 'varPhase';
    items_per(2).idx = [1 2];
    items_per(2).target = phase12_target;

    % Set control parameters
    items_controlled = [4 5];
    % Set stepsize
    % max_stepsize = abs([(period_target-StateDerived.period)/(2*pi) phase12_target-StateDerived.varPhase(1,2)])/200;
    max_stepsize = [0.1 0.1];
    errBound = 1e-6;

    Modtask = FMAM_ODE(sys,[],baseStateView,items_per,items_controlled,max_stepsize,errBound, ...
        'derivatives', derivatives);
    Modtask.psiUpdateMode = true;
    Modtask.refreshPsiModeReferences();
    Modtask.needLog = true;

    tic
    Modtask.fit()
    Modtask.step()
    elapsedTime = toc;  %
    disp(['Computing time: ', num2str(elapsedTime), ' seconds']);
    lambdaReached = Modtask.continuationStatus.lambda;
    if abs(Modtask.continuationStatus.lambda - 1) > 1e-12
        warning('Modulation_phase:ContinuationIncomplete', ...
            'Phase continuation stopped at lambda=%.6g for phase target %.6g.', ...
            Modtask.continuationStatus.lambda, phase12_target);
    end

    StateView = Modtask.exportSolverView();
    StateDerived = Modtask.exportDerivedView();

    Parameters = reshape(StateView.params, 1, []);
    y0 = StateDerived.TS_var(1, :).';
    searchWindow = max(10 * StateDerived.period, StateDerived.t(end));
    orbitOptions.tspan = [0, searchWindow];
    poResult = extract_periodic_orbit(odeFunc, y0, Parameters, orbitOptions);
    if ~poResult.has_orbit
        error('Modulation_phase:PeriodicOrbitGenerationFailed', ...
            'Periodic-orbit extraction failed for phase target %.6g (%s).', ...
            phase12_target, poResult.message);
    end

    nexttile
    hold on
    grid on
    box on
    plot(poResult.orbit_t, poResult.orbit_y(:,1), 'LineWidth', 1.2)
    plot(poResult.orbit_t, poResult.orbit_y(:,2), 'LineWidth', 1.2)
    plot(poResult.orbit_t, poResult.orbit_y(:,3), 'LineWidth', 1.2)
    title(sprintf('phase12 = %g, \\lambda = %.3g', phase12_target, lambdaReached))
    xlabel('t')
    ylabel('State')
    legend({'V_1', 'V_2', 'V_3'}, 'Location', 'best')
end
