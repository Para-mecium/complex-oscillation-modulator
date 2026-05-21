# complex-oscillation-modulator

## Reproducibility documentation

- `REPRODUCIBILITY.md`: repository-level guide for reproducing the numerical
  workflows from archived data or regenerated results.
- `FIGURE_WORKFLOW_MAP.md`: map from figure groups to scripts, input `.mat`
  files, and expected outputs.
- `ENVIRONMENT_AND_SEEDS.md`: pinned MATLAB environment, required toolboxes,
  MatCont note, and random-seed policy.

## Dependency
- `FMAM_ODE` itself only needs MATLAB.
- If you want to generate analytic Jacobians automatically, use `build_symbolic_derivatives.m`, which requires Symbolic Toolbox.

## Syntax
```matlab
derivatives = build_symbolic_derivatives(system, observables, numel(params));
continuationOptions = struct('predictorMode', 'auto', 'initialLambdaStep', 0.1);
stat = state(observables, params, t, TS_var, M, PV);
solverView = fmam_state_ops.solverViewFromState(stat);
task = FMAM_ODE(system, observables, solverView, items_per, items_controlled, ...
    [], errBound, 'derivatives', derivatives, ...
    'continuationOptions', continuationOptions);
```

## Description
This project is built to solve the complex oscillation modulation problem raised in article *Modulatability of complex oscillators*. The class `FMAM_ODE` computes the modulation path along which the properties of the periodic solution move towards the target value.

The derivative preparation step is now separated from the solver. `FMAM_ODE` expects a precomputed `derivatives` struct at construction time, so symbolic differentiation can be done once outside the class and reused across experiments.

## Input arguments
**system**: Functions of the ordinary differential equation, specified as a $1\times N$ cell where $N$ is the dimension of the ODE. Each components of 'system' are expected to be function handles.

**observables**: Functions that maps the state variables of ODE to the observables, specified as a $1\times n$ cell where $n$ is the number of observables. Each components of 'observables' are expected to be function handles.

**solverView**: Canonical solver-state struct storing the Fourier representation used by `FMAM_ODE`. If your workflow still starts from a `state` object, convert it with `fmam_state_ops.solverViewFromState(stat)` before constructing the task.

**items_per**: Struct array of modulation targets. Each entry must contain `prop`, `idx`, and `target`, for example `struct('prop','varAmp','idx',1,'target',10)`.

**items_controlled**: Row vector of parameter indices that are allowed to change during modulation.

**maxstepsize**: Scalar or row vector that caps each accepted continuation step in target space. Pass `[]` to remove this cap and let `continuationOptions` control the outer continuation step size directly.

**errBound**: Scalar or row vector used to detect already-satisfied targets and zero-length path components.

**continuationOptions**: Optional struct controlling the adaptive outer continuation and predictor-corrector policy. Supported fields include:
- `initialLambdaStep`, `initialSteps`, `minLambdaStep`, `maxLambdaStep`
- `growthFactor`, `shrinkFactor`, `backtrackingFactor`
- `easyIterations`, `slowIterations`, `maxFailures`
- `predictorMode` (`'auto'|'secant'|'quadratic'|'hermite'|'constant'`)
- `quadraticCurvatureThreshold`, `quadraticStepRatioBounds`
- `hermiteMaxExtrapolationRatio`, `predictorStepGrowthLimit`

**derivatives**: Precomputed analytic derivative cache passed as `'derivatives', derivatives`. It must provide:
- `derivatives.var(i,j).function`: Jacobians of the ODE right-hand side with respect to state variables and parameters.
- `derivatives.obs(i,j).function`: First derivatives of observables with respect to state variables.
- `derivatives.obs2(i,j,k).function`: Second derivatives of observables with respect to state variables.

Use `build_symbolic_derivatives` if you want to generate this cache from symbolic expressions before constructing `FMAM_ODE`.

Legacy examples under `Normal form/` still use older solver entrypoints and are not the recommended API examples for the current `FMAM_ODE` constructor.
