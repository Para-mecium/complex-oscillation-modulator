# complex-oscillation-modulator

MATLAB implementation for FMAM modulation examples, network experiments,
parameter-inference workflows, and figure-data generation associated with the
manuscript *Modulatability of complex oscillators* (under revision).

The class `FMAM_ODE` computes continuation paths in the parameter space along which properties of a
periodic solution move toward prescribed targets. Derivative preparation is
separated from the solver: symbolic differentiation is done once outside the
class, and the resulting cache is reused across experiments.

## Quick start

A typical modulation workflow proceeds in five stages. For a complete scripted
example, see `flexible_modulators/modulation_paths.m`.

### 1. Periodic orbit

Start from a stable periodic orbit of the ODE at the initial parameter values.
Most workflows obtain `(t, TS_var)` either by forward simulation with an event
detector, or through `PO_extract/extract_periodic_orbit.m`. The orbit must cover at least one full period after any
transient has been discarded.

### 2. Phase representation

Choose a primary variable `PV` (a state component or observable) and a Fourier
truncation order `M`. Build a `state` object from the sampled orbit, then
convert it to the canonical solver struct:

```matlab
stat = state(observables, params, t, TS_var, M, PV);
stat.updatePeriod();
stat.updateVar2();   % and updateObs2() when observables are used
solverView = fmam_state_ops.solverViewFromState(stat);
```

`solverView` stores the Fourier representation (`p_Psi`, `q_Psi`, `p_var`,
`q_var`, …) that `FMAM_ODE` operates on.

### 3. Derivative cache

Prepare analytic Jacobians before constructing the task. With the Symbolic
Toolbox:

```matlab
derivatives = build_symbolic_derivatives(system, observables, numel(params));
```

Alternatively, supply an equivalent struct of function handles with fields
`var`, `obs`, and `obs2` (see **Input arguments** below).

### 4. Modulation setup and solve

Define one modulation target per struct entry in `items_per`, and pair each
target with one controlled parameter index in `items_controlled`:

```matlab
items_per(1).prop = 'varAmp';
items_per(1).idx = 1;
items_per(1).target = 10;
items_per(2).prop = 'p_Psi';
items_per(2).idx = 1;
items_per(2).target = targetPeriod / (2 * pi);

items_controlled = [1, 2];   % one parameter index per target
accuracy = 1e-6;
continuationOptions = struct('predictorMode', 'auto', 'initialLambdaStep', 0.1);

task = FMAM_ODE(system, observables, solverView, items_per, items_controlled, ...
    [], accuracy, 'derivatives', derivatives, ...
    'continuationOptions', continuationOptions);

task.fit();    % Newton solve at the current target values
task.step();   % adaptive continuation toward the modulation targets

task.errBound = 1e-12;   % optional: tighten Newton tolerance after continuation
task.fit();
```

Supported values for `items_per(i).prop` are listed under **Modulation
targets** below.

### 5. Export and post-processing

Read results from the task object rather than rebuilding from raw coefficients
when possible:

```matlab
solverView = task.exportSolverView();
derivedView = task.exportDerivedView();   % period, amplitudes, extrema, ...
```

Use `task.continuationStatus` and `task.logs` (when `needLog` is enabled) to
inspect the continuation trace.

## Reproducing manuscript figures

Each figure group in the repository follows three stages:

1. Plotting-data generation with the corresponding `build_*`, `run_*`, or
   inference scripts.
2. FMAM modulation, parameter inference, continuation, or repeated simulation.
3. Plotting from the corresponding `.mat` outputs.

Plotting `.mat` files are generated outputs and are not expected to be tracked
by Git. Before running a `draw_*` script, generate its inputs with the production
scripts listed in `FIGURE_WORKFLOW_MAP.md` and the relevant `WORKFLOW.md` file.
Some full workflows contain computationally expensive sweeps; their documented
software requirements and random-seed policies apply during regeneration.

| Goal | Where to look |
|---|---|
| Map figure groups to scripts, inputs, and outputs | `FIGURE_WORKFLOW_MAP.md` |
| Pinned MATLAB/toolbox versions and random-seed policy | `ENVIRONMENT_AND_SEEDS.md` |

To reproduce a panel, run the listed `build_*`, `run_*`, or inference scripts
first, then run the corresponding `draw_*` script.

## Repository layout

| Path | Role |
|---|---|
| `FMAM_ODE.m`, `state.m`, `fmam_state_ops.m` | Core modulation solver and periodic-orbit representation |
| `build_symbolic_derivatives.m` | Optional symbolic derivative cache builder |
| `solve_generic_newton.m`, `solve_regularized_linear_system.m` | Newton and linear solver support |
| `PO_extract/` | Periodic-orbit extraction  |
| `flexible_modulators/`, `Circadian/`, `Longevity/`, `RLT_circuit/`, `network_modulatability/`, `Normal form/`, `complexity_analyses/` | Manuscript figure and supplementary workflows |

Add the repository root (and `PO_extract/` when needed) to the MATLAB path before
running scripts in any subdirectory.

## Dependency

- `FMAM_ODE` itself only needs MATLAB.
- At construction time, `FMAM_ODE` requires a precomputed `derivatives` struct.
  Generate it with `build_symbolic_derivatives.m` (Symbolic Math Toolbox), or
  supply an equivalent struct of function handles yourself.
- Specific figure workflows may additionally require Parallel Computing Toolbox,
  Statistics and Machine Learning Toolbox, or Signal Processing Toolbox; see
  `ENVIRONMENT_AND_SEEDS.md`.
- Some periodic-orbit workflows in `PO_extract/` can use [MatCont](https://sourceforge.net/projects/matcont/) as an optional backend. MatCont is not bundled in this repository; place it under `MatCont7p6/` at the repository root when needed.

## Modulation targets

Each entry in `items_per` must contain `prop`, `idx`, and `target`.

Supported `prop` values:

| `prop` | Quantity modulated |
|---|---|
| `params` | Model parameter (direct target; index must also appear in `items_controlled`) |
| `p_Psi`, `q_Psi` | Fourier coefficients of `Psi(phi) = dt/dphi` |
| `p_var`, `q_var` | Fourier coefficients of a state variable |
| `varPhiMax`, `varPhiMin` | Phase of a state-variable extremum |
| `obsPhiMax`, `obsPhiMin` | Phase of an observable extremum |
| `varAmp`, `obsAmp` | Amplitude of a state variable or observable |
| `varMax`, `varMin`, `obsMax`, `obsMin` | Extremum value |
| `varPhase` | Phase difference between two state variables or observables |

`idx` format depends on `prop`: Fourier-coefficient targets use `[i, j]` (variable/observable index and harmonic index); `varPhase` uses `[i, j]` for the two components whose phase difference is targeted; scalar properties such as `varAmp` use the component index alone. See the header comment in `FMAM_ODE.m` for details.

## Input arguments

**system**: Functions of the ordinary differential equation, specified as a $1\times N$ cell where $N$ is the dimension of the ODE. Each components of 'system' are expected to be function handles.

**observables**: Functions that maps the state variables of ODE to the observables, specified as a $1\times n$ cell where $n$ is the number of observables. Each components of 'observables' are expected to be function handles.

**solverView**: Canonical solver-state struct storing the Fourier representation used by `FMAM_ODE`. If your workflow still starts from a `state` object, convert it with `fmam_state_ops.solverViewFromState(stat)` before constructing the task.

**items_per**: Struct array of modulation targets. Each entry must contain `prop`, `idx`, and `target`, for example `struct('prop','varAmp','idx',1,'target',10)`.

**items_controlled**: Row vector of parameter indices, one per modulation target. Its length must equal `numel(items_per)`, indices must be unique, and each entry must lie in `1:dimParams`.

**maxstepsize**: Scalar or row vector that caps each accepted continuation step in target space. Pass `[]` to remove this cap and let `continuationOptions` control the outer continuation step size directly.

**accuracy**: Scalar or row vector (7th positional argument) used to detect already-satisfied targets and zero-length path components during continuation. Targets whose remaining distance falls below `accuracy` are treated as inactive.

**errBound**: Optional name-value pair (default `1e-8`) controlling Newton convergence tolerance. This is separate from the 7th positional `accuracy` argument; scripts often use a looser `accuracy` during continuation and then set a tighter `task.errBound` before a final `fit()`.

**continuationOptions**: Optional struct controlling the adaptive outer continuation and predictor-corrector policy. Supported fields include:
- `initialLambdaStep`, `initialSteps`, `minLambdaStep`, `maxLambdaStep`
- `growthFactor`, `shrinkFactor`, `backtrackingFactor`
- `easyIterations`, `slowIterations`, `maxFailures`
- `predictorMode` (`'auto'|'secant'|'quadratic'|'hermite'|'constant'`)
- `quadraticCurvatureThreshold`, `quadraticStepRatioBounds`
- `hermiteMaxExtrapolationRatio`, `predictorStepGrowthLimit`
- `conditionStopEnabled`, `conditionStopRcond`

**derivatives**: Precomputed analytic derivative cache passed as `'derivatives', derivatives`. It must provide:
- `derivatives.var(i,j).function`: Jacobians of the ODE right-hand side with respect to state variables and parameters.
- `derivatives.obs(i,j).function`: First derivatives of observables with respect to state variables.
- `derivatives.obs2(i,j,k).function`: Second derivatives of observables with respect to state variables.

Use `build_symbolic_derivatives` if you want to generate this cache from symbolic expressions before constructing `FMAM_ODE`.

## License

The software in this repository is licensed under the [MIT License](LICENSE).

## Citation

<!-- TODO: add bibliographic details once the manuscript is accepted -->

If you use this code, please cite:

> *[Modulatability of complex oscillators]* — citation details to be added after publication.
