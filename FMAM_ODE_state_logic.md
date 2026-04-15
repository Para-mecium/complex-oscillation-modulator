# Logic Extraction for `FMAM_ODE.m` and `state.m`

This note reflects the current implementation of the two core MATLAB classes in `FMAM_code`:

- `state.m`: preprocesses one sampled periodic orbit into the FMAM phase representation and stores both Fourier coefficients and derived orbit properties.
- `FMAM_ODE.m`: solves modulation/continuation problems directly on those Fourier coefficients, auxiliary extremum phases, and selected model parameters.

The code is not a direct ODE integrator. It assumes one periodic orbit is already known as sampled data, then rewrites that orbit in a phase variable `phi` so that one chosen primary variable is close to a cosine. Modulation is then performed in coefficient space.

## 1. High-level workflow

The intended workflow is:

1. Start from an already computed periodic orbit `t, TS_var`.
2. Choose a primary variable `PV`.
   - `PV.name = 'var'` uses one state component as the primary variable.
   - `PV.name = 'obs'` uses one observable as the primary variable.
3. Build a `state` object.
   - normalize the input time axis;
   - remove a duplicated periodic endpoint if present;
   - build Fourier coefficients for `Psi(phi) = dt/dphi`, state variables, and observables;
   - compute derived quantities such as period, extrema, amplitudes, and phase differences.
4. Build an `FMAM_ODE` object.
   - define which quantities should be modulated in `items_perturb`;
   - define which model parameters are allowed to move in `items_controlled`;
   - precompute symbolic derivatives of the ODE and observables.
5. Call `fit()` to solve the inner linearized correction problem for the current targets.
6. Call `step()` to advance the targets toward their final requested values by outer continuation.
7. Prefer `task.exportSolverView()` and `task.exportDerivedView()` for downstream post-processing.
8. Rebuild a rich `state` object only for compatibility paths that still need object-style access.

The canonical runtime boundary is now `solverView + derivedView`. The `state` class remains a compatibility layer around those views rather than the source of truth for solver transactions.

## 2. Representation used by `state.m`

### 2.1 Core idea

The periodic orbit is represented in a transformed phase variable `phi in [0, 2*pi)` rather than directly in physical time `t`.

The main unknowns are:

- `p_Psi`, `q_Psi`: Fourier coefficients of `Psi(phi) = dt/dphi`
- `p_var`, `q_var`: Fourier coefficients of the state variables
- `p_obs`, `q_obs`: Fourier coefficients of the observables

With truncation order `M`, the reconstruction is:

```text
Psi(phi) = sum_{k=0}^{M-1} p_Psi(k) cos(k phi) + sum_{k=1}^{M-1} q_Psi(k) sin(k phi)

x_i(phi) = sum_{k=0}^M p_var(k,i) cos(k phi) + sum_{k=1}^M q_var(k,i) sin(k phi)
```

The stored period is

```text
T = 2*pi*p_Psi(1)
```

because the constant Fourier mode of `Psi` is the average value of `dt/dphi`.

### 2.2 Constructor logic

`state(obs, Params, t, TS_var, M, PV)` currently does the following:

1. Normalize the observable container:
   - `[]` becomes `{}`;
   - a cell input is reshaped to a row cell array.
2. Normalize `t` and `Params` to column/row vectors.
3. Call `normalizePeriodicInputs(...)`.
   - shift the input time axis so it starts at zero;
   - if the sampled orbit includes a repeated periodic endpoint, drop the last sample.
4. Validate the inputs with `validateInputs(...)`.
5. Evaluate observables on the sampled orbit with `getObs(...)`.
6. Call `fourierCoeffs(...)` to build the FMAM phase representation.
7. Store:
   - the Fourier coefficients;
   - the primary-variable cosine parameters `a` and `b`;
   - the primary variable spec `PV`.
8. Call `refreshDerivedState()`, which currently:
   - checks the time-map invariant;
   - updates equal-time Fourier coefficients `p_var_origin`, `q_var_origin`;
   - updates variable extrema/amplitudes/phases;
   - updates observable extrema/amplitudes/phases;
   - updates the period.

### 2.3 Input normalization and validation

Compared with the older code path, `state` now enforces a clearer input contract:

- `t` may start at any origin; the constructor internally shifts it to zero.
- a repeated periodic endpoint is accepted and trimmed automatically.
- `t` must be strictly increasing after trimming.
- `TS_var` must match the length of `t` and contain only finite values.
- the period span `t(end)` after normalization must be positive.
- `PV` must name a valid state or observable index.

### 2.4 What `fourierCoeffs(...)` does

This is the key preprocessing step.

1. Normalize the input time vector again to `t - t(1)`.
2. Define the period as `T = t(end)`, which is now the time span of one period after normalization.
3. Resample the orbit and observables onto a uniform time grid.
4. Extract the chosen primary variable `X`.
5. Rotate the sampled orbit so the primary variable starts at its maximum.
6. Estimate
   - `a = (max(X)-min(X))/2`
   - `b = (max(X)+min(X))/2`
7. Split the primary waveform into two usable branches:
   - max-to-min branch `1:I2`
   - min-to-max branch `I2:end`
8. Match `a*cos(phi)+b` against those two branches to build a phase-to-time map `t(phi)`.
9. Differentiate `t(phi)` with a finite-difference matrix `Atrans` to obtain sampled `Psi(phi)`.
10. Interpolate the state and observable samples onto the induced `phi` grid.
11. FFT-project `Psi` to order `M-1`, and states/observables to order `M`.
12. Enforce the primary-variable normalization:
   - if `PV` is a state variable, its Fourier series is forced to `b + a cos(phi)`;
   - if `PV` is an observable, the chosen observable is forced to have only the first cosine mode.
13. Validate the time map:
   - sampled `Psi(phi)` must stay finite and strictly positive on the sampling grid;
   - exact interval increments from integrating `Psi` over each phase cell must also stay strictly positive.

If the primary waveform does not admit a usable max-to-min and min-to-max branch, the code raises `state:InvalidPrimaryWaveform`.

### 2.5 Dependent quantities

`state` reconstructs several arrays on demand:

- `phi`: a uniform phase grid of length `LphiConst`
- `Psi`: the reconstructed trigonometric polynomial for `dt/dphi`
- `TS_var`: the reconstructed state trajectory on the uniform `phi` grid
- `TS_obs`: the reconstructed observable trajectory on the same grid
- `t`: a reconstructed time grid obtained by cumulative summation of exact phase-cell integrals of `Psi`

Important detail:

- `Atrans` still exists and is used in preprocessing to estimate sampled `Psi` from `t(phi)`.
- but the runtime `t` property is no longer obtained by solving an inverse finite-difference system.
- instead, `get.t` first checks sampled `Psi > 0`, then builds time increments from exact trigonometric integrals with `Trintegration(...)`.

This avoids the old inconsistency where a finite-difference reconstruction of `t` could become non-monotone even though the underlying integrated `Psi` remained positive.

### 2.6 Time-map invariants

`state` now exposes three related pieces of logic:

- `assertPositivePsi(...)`: sampled `Psi(phi)` must remain finite and strictly positive.
- `assertPositiveTimeIncrements(...)`: exact interval integrals of `Psi` over the phase grid must remain finite and strictly positive.
- `assertTimeMapInvariant()`: runs both checks on the current stored coefficients.

These invariants are enforced:

- during preprocessing in `fourierCoeffs(...)`;
- during `refreshDerivedState()`;
- when `t` is requested.

They are not currently enforced inside the inner Newton loop of `FMAM_ODE.fit()`.

### 2.7 Derived orbit properties

`updateVar2()` and `updateObs2()` compute:

- maximum phases: `varPhiMax`, `obsPhiMax`
- minimum phases: `varPhiMin`, `obsPhiMin`
- amplitudes: `varAmp`, `obsAmp`
- extrema values: `varMax`, `varMin`, `obsMax`, `obsMin`
- phase-difference matrices: `varPhase`, `obsPhase`

The extrema are found by `FindExtreme(...)`, which:

1. samples the trigonometric polynomial on a grid;
2. locates discrete maxima and minima;
3. checks the derivative residual at those locations;
4. doubles the sampling grid until the residual is small enough or the iteration cap is hit.

Phase differences are computed in physical time by integrating `Psi` between extremum phases with `Trintegration(...)`.

### 2.8 `p_var_origin` and `q_var_origin`

`updateFCOrigin()` computes a separate equal-time Fourier description of the orbit:

- it uses the current `t` grid and reconstructed `TS_var`;
- it resamples each state variable onto an equal-time grid;
- it stores the resulting coefficients in `p_var_origin` and `q_var_origin`.

These coefficients are not the main FMAM unknowns. They are an auxiliary time-domain Fourier description derived from the phase-domain solution.

### 2.9 Snapshot-based rebuild support

`state` now provides:

- `state.fromSolverSnapshot(obs, snapshot)`

The snapshot includes the coefficient blocks and derived scalar/matrix properties needed for compatibility rebuilds. `FMAM_ODE` now keeps its runtime continuation state directly as `solverView` structs, and rich `state` objects are rebuilt explicitly when needed.

## 3. Representation used by `FMAM_ODE.m`

### 3.1 Main responsibility

`FMAM_ODE` treats modulation as a continuation problem:

- slightly move the current target values toward the final requested targets;
- solve a linearized correction problem for the orbit and selected parameters;
- accept or reject that continuation step;
- repeat until all targets are reached.

### 3.2 Important inputs

- `sys`: cell array of ODE right-hand-side functions
- `obs`: cell array of observable functions
- `initialSolverView`: canonical `solverView` struct storing the runtime continuation state
- `items_perturb`: list of quantities to modulate
- `items_controlled`: indices of parameters allowed to move
- `maxstepsize`: cap for each outer continuation step
- `accuracy`: tolerance used by `autostepsize()`

Each entry of `items_perturb` has:

- `prop`: property name to modulate
- `idx`: component index
- `target`: final desired value

The constructor now validates:

- every target property name against the supported target set;
- target indices against the actual dimensions of the current `state`;
- `items_controlled` for uniqueness and parameter-range validity;
- `numel(items_controlled) == numel(items_perturb)`.

### 3.3 Quantities that can be targeted

The current solver branches support these target types:

- direct unknowns:
  - `params`
  - `p_Psi`, `q_Psi`
  - `p_var`, `q_var`
  - `varPhiMax`, `varPhiMin`
  - `obsPhiMax`, `obsPhiMin`
- derived targets:
  - `varAmp`, `varMax`, `varMin`
  - `obsAmp`, `obsMax`, `obsMin`
  - `varPhase`

`obsPhase` is still exposed at the `state` level but is not implemented as a modulation target in the current `FMAM_ODE` target logic.

### 3.4 Constructor logic

The constructor:

1. normalizes the `system` and `observables` containers;
2. applies optional name-value overrides;
3. stores the system, observables, current state, and solver tolerances;
4. derives dimensions from the `state` object;
5. validates the modulation setup;
6. initializes:
   - `targetCurr` from the current orbit, via `initialTargetValues()`;
   - `p_Psi_init`, `q_Psi_init` when `isPsiUpdated == false`;
7. uses Symbolic Toolbox to build:
   - Jacobians of every ODE component with respect to state variables and parameters;
   - first derivatives of observables;
   - second derivatives of observables.

## 4. Outer continuation logic: `step()`

### 4.1 `autostepsize()`

`autostepsize()` is now only a helper for assembly-oriented tests.

It uses the same continuation path as `step()`:

1. builds a scalar `lambda` path from current targets to final targets;
2. computes an initial `dlambda`;
3. sets `stepsize = target(lambda + dlambda) - targetCurr`.

The production continuation loop is implemented in `step()`.

### 4.2 `items_per_curr`

`items_per_curr` is a dependent property that simply returns `targetCurr`.

This matters because the outer continuation now distinguishes:

- the final requested values in `items_perturb(i).target`;
- the current accepted continuation targets in `targetCurr`.

### 4.3 Current transactional `step()` loop

The outer loop is an accepted-state-driven predictor-corrector continuation:

1. build continuation path `target(lambda)`, with `lambda in [0,1]`;
2. initialize from an accepted snapshot at `lambda = 0`;
3. propose `lambdaTrial = lambda + dlambdaTrial`;
4. set tentative `targetCurr = target(lambdaTrial)` and prepare predictor state (`constant/secant/quadratic/hermite`);
5. run `fit()` as the corrector;
6. if correction fails, rollback to the latest accepted snapshot, shrink `dlambda`, and retry;
7. if correction succeeds, commit the state and update `dlambda` using Newton iteration feedback;
8. stop when `lambda >= 1`, or when `dlambda` falls below `minLambdaStep`, or failures exceed `maxFailures`.

Only accepted continuation points are logged.

## 5. Inner correction logic: `fit()` and `res()`

### 5.1 `res()`

`res()` returns four residual magnitudes:

1. `res_sys`: ODE residual in Fourier space
2. `res_var_phi`: residual of the variable-extremum conditions
3. `res_obs_phi`: residual of the observable-extremum conditions
4. `res_target`: mismatch between the current orbit and the current continuation targets `targetCurr`

The ODE residual is built from

```text
- d x / dphi + Psi(phi) * f(x(phi), params)
```

sampled on a dense phase grid and then transformed by FFT.

The target residual depends on the target type:

- `varAmp` uses max minus min
- `varMax`, `varMin` evaluate the variable at variable extremum phases
- `obsAmp`, `obsMax`, `obsMin` evaluate the observable at observable extremum phases
- `varPhase` integrates `Psi` between variable maxima

### 5.2 `fit()`

`fit()` is the inner corrector entry for the current `targetCurr`.

It now returns a struct:

- `result.converged`
- `result.iterations`
- `result.finalError`

`fit()` delegates the iteration to `solve_generic_newton` through `runNewtonSolve`.
Termination, line search/regularization, and diagnostics all come from the generic Newton solver and `newtonOptions`.

Important current behavior:

- `fit()` no longer performs `Psi`-feasibility or time-map validation itself.
- it only monitors the algebraic residuals returned by `res()`.
- time-map validation remains in `state`, and is triggered when the code explicitly rebuilds derived quantities or requests the physical-time grid.

## 6. Linearized solve in `oneIter()`

### 6.1 Unknown vector

The solver forms one linear system

```text
A * increments = res
```

with unknown blocks:

- `params`
- `p_Psi`, `q_Psi`
- `p_var`, `q_var`
- `varPhiMax`, `varPhiMin`
- `obsPhiMax`, `obsPhiMin`

An index map is built so each property can be written into the correct block of the global vector.

### 6.2 Equation blocks

The rows of `A` come from several blocks.

#### 6.2.1 Block A: Incremental ODE equations

For each state equation:

1. linearize `-x_phi + Psi*f(x,p)`;
2. transform the coefficient functions with FFT;
3. keep cosine and sine modes up to order `M`.

This couples parameter increments, `Psi` increments, and state Fourier increments.

#### 6.2.2 Block B: Extremum constraints

For each variable/observable maximum or minimum needed by the current targets:

- enforce zero derivative at that extremum phase;
- include derivatives with respect to Fourier coefficients and the extremum phase itself.

That is why `varPhiMax`, `varPhiMin`, `obsPhiMax`, and `obsPhiMin` are part of the unknown vector.

#### 6.2.3 Block C: Fixed-parameter constraints

Any parameter not listed in `items_controlled` receives a zero-increment constraint.

So only the selected parameters are allowed to move.

#### 6.2.4 Block D: Target equations

For each entry in `items_perturb`, the solver adds one linearized equation for the current continuation target `targetCurr(j)`:

- direct coefficient targets use identity rows;
- amplitude/max/min/phase targets use the appropriate analytic linearization.

#### 6.2.5 Block E: Gauge / normalization constraints

The code removes the phase/gauge freedom in one of two ways:

- `isPsiUpdated == false`:
  higher Fourier modes of `Psi` are constrained to stay at their initial values.
- `isPsiUpdated == true`:
  the primary-variable normalization is enforced instead.
  - if `PV` is a state variable, higher harmonics and sine modes of that state are forced to zero;
  - if `PV` is an observable, the corresponding observable normalization is enforced through observable derivatives.

### 6.3 Solve path

The linearized system is assembled in `assemble_newton_linear_system`, then solved and applied through `solve_generic_newton` callbacks:

1. `linearize()` builds `A`, residual, and unknown index map;
2. `applyIncrement()` updates unknown blocks and validates the time-map invariant;
3. `measure()` returns residual diagnostics and convergence status.

Any rejected Newton trial is rolled back by the snapshot/restore contract inside the generic solver.

## 7. Static helper methods in `FMAM_ODE`

The static methods implement the analytic pieces used by the linearization:

- `Vec_CS`: builds cosine/sine basis matrices
- `Trintegration`: integrates a trigonometric polynomial between two phases
- `delta_coe_system`: derivatives of the ODE residual with respect to parameters and Fourier coefficients
- `residue_system`: raw ODE residual
- `delta_coe_phi_var`, `residue_phi_var`: variable-extremum condition and its linearization
- `delta_coe_obs_phi`, `residue_phi_obs`: observable-extremum condition and its linearization
- `delta_coe_var_target`, `var_target_curr`: variable max/min/amplitude target formulas
- `delta_coe_obs_target`, `obs_target_curr`: observable target formulas
- `delta_coe_state_phase`: linearization of time-phase difference between two maxima

These methods provide the closed-form pieces needed by `oneIter()` without re-deriving each target equation manually.

## 8. Practical interpretation

The current solver can be read as:

1. encode one periodic orbit in a phase coordinate where one primary quantity is normalized;
2. treat that orbit as a finite-dimensional vector of Fourier coefficients plus auxiliary extremum phases;
3. store a current accepted continuation target vector `targetCurr`;
4. move `targetCurr` one outer step toward the user-requested final targets;
5. solve a damped linearized correction problem for coefficients, extrema, and selected parameters;
6. accept or roll back the entire continuation step;
7. repeat until all requested targets are reached.

This is a continuation method on orbit properties, not a time integration method.

## 9. Important implementation notes

- The code assumes an initial periodic orbit is already known.
- Symbolic Toolbox is required because the `FMAM_ODE` constructor differentiates the ODE and observable functions automatically.
- `state` now accepts nonzero time origins and repeated periodic endpoints.
- `state` enforces a time-map invariant when it preprocesses an orbit or reconstructs the physical-time grid.
- `FMAM_ODE.fit()` does not currently enforce that invariant during the Newton loop.
- `fit()` returns a result struct and does not throw on ordinary non-convergence.
- `step()` is transactional and rolls back both the stored orbit state and `targetCurr` when a continuation step fails.
- logging is append-only for accepted continuation steps via `appendAcceptedLog()`.
- `items_perturb` can target both raw Fourier quantities and higher-level orbit features such as amplitude, extrema, and phase lag.
- `items_controlled` only refers to model parameters. Other orbit quantities are adjusted indirectly through the linearized solve.

## 10. Short mental model

`state.m` answers:

> How do we encode one sampled periodic orbit as FMAM Fourier data in a normalized phase coordinate, while keeping a consistent physical-time map?

`FMAM_ODE.m` answers:

> If I slightly change the desired orbit properties, how do I update parameters, coefficients, and extremum phases so the orbit approximately remains a periodic solution, and how do I roll back safely when one continuation step fails?
