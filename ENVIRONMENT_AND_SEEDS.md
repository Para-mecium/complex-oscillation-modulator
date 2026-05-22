# Environment And Seed Policy

This file pins the software environment used for the current repository checkout
and documents how stochastic or randomized workflows handle random seeds.

## MATLAB Environment

The repository was checked with:

```text
MATLAB version: 24.2.0.2806996 (R2024b) Update 3
Operating system: macOS 26.4.1
Java: 11.0.25+9-LTS
```

Installed MATLAB products available in the checked environment:

```text
MATLAB R2024b
MATLAB Coder R2024b
MATLAB Compiler R2024b
MATLAB Compiler SDK R2024b
MATLAB Report Generator R2024b
Parallel Computing Toolbox R2024b
Signal Processing Toolbox R2024b
Statistics and Machine Learning Toolbox R2024b
Symbolic Math Toolbox R2024b
```

The core `FMAM_ODE` solver only requires MATLAB. Specific workflows use the
following additional components:

- Symbolic Math Toolbox: `build_symbolic_derivatives.m` and scripts that
  regenerate analytic derivative caches.
- Parallel Computing Toolbox: `parfor` workflows, including repeated SDE
  simulations and large baseline/parameter sweeps.
- Statistics and Machine Learning Toolbox: sampling, distribution summaries, and
  statistical helper functions used by selected plotting or analysis scripts.
- Signal Processing Toolbox: PSD-related analysis in stochastic circadian
  workflows.

## Seed Policy

All stochastic or randomized workflows should set `rng(seed, 'twister')` through
the script-level seed value or through a `seed` field stored in the workflow
options. Scripts that run repeated simulations store the seed vectors in the
corresponding `.mat` output.

The deterministic FMAM continuation scripts do not use random seeds unless the
particular model setup samples a network, initial condition, noisy observation,
or stochastic trajectory.

## Documented Seed Locations

| Workflow | Seed handling |
|---|---|
| Longevity SDE representative panel | `Longevity/build_etd_fig2c_data.m` uses `seed = 1`; the saved output stores `seed` and `Options`. |
| Longevity SDE ensemble panel | `Longevity/build_etd_fig2d_data.m` uses `distributionSeeds = 1001:1100`; the saved output stores `distributionSeeds` and `distributionOptions`. |
| Longevity SDE simulator | `Longevity/SDE_simulation/Longevity_SDE.m` defaults to `opts.seed = 1` and calls `rng(double(opts.seed), 'twister')`. |
| Circadian SDE representative panel | `Circadian/build_fig5d_sde_representative_data.m` uses `representativeSeeds = [0, 0, 0]`; the saved output stores `representativeSeeds`. |
| Circadian repeated SDE statistics | `Circadian/build_fig5d_sde_repeat_data.m` uses `repeatSeeds = 1:100`; the saved output stores `repeatSeeds`. |
| Circadian SDE simulator | `Circadian/SDE_simulation/Circadian_SDE.m` defaults to `opts.seed = 1` and calls `rng(double(opts.seed), 'twister')`. |
| Network generation | `network_modulatability/+networkexp/build_network_matrix.m` temporarily sets `rng(seed, 'twister')` and restores the previous RNG state when finished. |
| Network single/ergodic experiments | See [Network modulatability seeds](#network-modulatability-seeds) below. |
| Circuit sensitivity to initial data | See [Circuit sensitivity: initial data](#circuit-sensitivity-initial-data) below. |
| Circuit sensitivity to measurement noise | See [Circuit sensitivity: measurement noise](#circuit-sensitivity-measurement-noise) below. |
| Baseline comparison | `RLT_circuit/parameter_inference/baseline_comparison/run_baseline_comparison.m` uses `randomSeeds = 1:100`; individual methods call `rng(seed, 'twister')`. |
| Plot jitter only | Scripts such as `draw_baseline_comparison.m` and `draw_sensitivity_to_init_data_violin.m` use fixed `rng(1, 'twister')` before visual scatter jitter; this affects marker placement but not numerical results. |

### Network modulatability seeds

All network workflows share one **master seed** in the top-level script: `settings.randomSeed` (default `1` in `FHN_network.m`, `GRN_network.m`, `FHN_network_ergodic.m`, and `GRN_network_ergodic.m`). Child seeds are **not** drawn from MATLAB’s global RNG; they are computed deterministically by `network_modulatability/+networkexp/derive_seed.m`, which mixes `masterSeed` with integer case indices (modulus `4294967291`, LCG-style update). Each derived seed is then passed explicitly to `rng(seed, 'twister')` inside `build_network_matrix.m`, `sample_perturbation.m`, or `sample_perturbation_from_sources.m` (with the previous RNG state restored afterward where applicable).

**Single experiment** (`networkexp.run_single_experiment`, launched from e.g. `FHN_network.m`):

| Role | Derivation | Used for |
|---|---|---|
| Baseline topology | `derive_seed(randomSeed, 1)` | `build_network_matrix` → baseline adjacency/weights |
| Perturbation edges | `derive_seed(randomSeed, 2, nPerturbedEdges)` | `sample_perturbation` (internally splits into `derive_seed(seed, 11)` for source nodes and `derive_seed(seed, 29)` for target nodes when modes are random) |

Outputs under `{Model}_{netName}/`: baseline file `TS_init_*.mat` (`baselineData.config`, `baselineData.network.metadata.seed`) and perturbation file `TS_per_{n}_*.mat` (`perturbedData.perturbation` with `seed`, `sourceSeed`, `targetSeed`; perturbed network metadata reuses the baseline network seed because the topology is the same matrix with a mask added).

**Ergodic sweep** (`networkexp.run_ergodic_experiment`, launched from e.g. `FHN_network_ergodic.m`):

| Role | Derivation | Used for |
|---|---|---|
| Baseline topology (per net type) | `derive_seed(randomSeed, idxNet, 1)` | One baseline network per entry in `settings.nets` |
| Source-node order (per sequence) | `derive_seed(randomSeed, idxNet, idxSourceSequence, 101)` | `randperm` when `perturbedNodeSelectionMode` is random; prefix mode uses `1:nPerturbedEdges` without RNG |
| Repeat sample ID | `derive_seed(randomSeed, idxNet, idxSourceSequence, nPerturbedEdges, attemptCount, 1)` | Stored as `sampleSeed` |
| Perturbation draw | `…, attemptCount, 2)` | `sample_perturbation_from_sources` → `perturbationSeed` |
| Network redraw (optional) | `…, attemptCount, 3)` | Only if `settings.repeatMode` is `'independentNetworkDraws'`; otherwise `networkSeed` equals the baseline seed and `baselineOmega` is reused |

For each combination of network type, source sequence, and `nPerturbedEdges = 1:settings.network.N`, the sweep runs `settings.numRepeat` attempts (or up to `numRepeat * maxAttemptsMultiplier` when `repeatAccountingMode` is `'untilSuccesses'`). Saved under `Ergodic data/N = {N}/` as one `.mat` per net (and per `source_seq` when not in prefix mode). Top-level fields: `experiment.config.randomSeed`, `experiment.seedPolicy` (`masterSeed`, `baselineNetworkSeed`, `sourceSequenceSeed`, `repeatMode`, `repeatAccountingMode`), and per-attempt records in `experiment.repeats(*).samples(*)` with `networkSeed`, `perturbationSeed`, `sampleSeed`, plus packed `perturbation` / `networkSummary` seeds.

To reproduce a case, set the same `settings.randomSeed` and experiment options in the entry script, rerun the corresponding `*_network.m` or `*_network_ergodic.m`, or regenerate one draw by calling `derive_seed` with the same index arguments and then `build_network_matrix` / `sample_perturbation*` with that child seed.

### Circuit sensitivity: initial data

Script: `RLT_circuit/parameter_inference/sensitivity_to_init_data/build_init_data_files.m`.

Each scale uses its own **base seed**; every sampling attempt uses `sampleSeed = baseSeed + attemptCount` with `rng(sampleSeed, 'twister')` before drawing four parameters (`sampledIdx = [1 4 5 6]`) uniformly in the scaled box. Only attempts that yield a periodic orbit via `circuit_forward_orbit` are saved; the variable `seed` written to disk equals that attempt’s `sampleSeed` (not the file index).

| `scale` | `baseSeed` | Attempt seeds tried (in order) | Output directory |
|---|---|---|---|
| `1` | `0` | `1:100`  | `init_data_files/scale_1/` |
| `0.5` | `1000` | `1001:1100`  | `init_data_files/scale_0p5/` |
| `0.75` | `2000` | `2001:2100`  | `init_data_files/scale_0p75/` |
| `1.25` | `3000` | `3001:3100`  | `init_data_files/scale_1p25/` |
| `1.5` | `4000` | `4001:4100`  | `init_data_files/scale_1p5/` |

Other fixed settings in the builder: `N = 100` accepted files per scale, `maxAttempts = 100`, parameters initialized from `parameter_inference/initData_ODE.mat` with non-sampled indices fixed to the ODE baseline. Accepted files are named `initData_001.mat` … `initData_100.mat`; the `seed` field inside each file is whichever `sampleSeed` produced that acceptance (typically the first 100 successful seeds in increasing `attemptCount`, but failed orbits skip indices). Downstream `build_sensitivity_to_init_data.m` currently processes scales `1.25` and `1.5` only; analysis/plot scripts use `[0.5 0.75 1 1.25 1.5]` where noted in their headers.

### Circuit sensitivity: measurement noise

Script: `RLT_circuit/parameter_inference/sensitivity_to_measurement_noise/build_noisy_init_data_files.m`.

Clean trajectory comes from `parameter_inference/initData_circuit.mat` (two periods concatenated). For each noise level, file `i = 1:100` uses the same seed schedule:

| `noiseLevel` | Per-file `seed` | Files |
|---|---|---|
| `0.01` | `1:100` | `noisy_init_data_files/noise_level_0p01/initData_001.mat` … `initData_100.mat` |
| `0.02` | `1:100` | `noisy_init_data_files/noise_level_0p02/initData_*.mat` |
| `0.05` | `1:100` | `noisy_init_data_files/noise_level_0p05/initData_*.mat` |

Noise is additive: `yNoisy = yTwoPeriod + noiseLevel * randn(size(yTwoPeriod))` immediately after `rng(seed, 'twister')`. Each `.mat` stores `seed`, `noiseLevel`, Lomb-estimated `period`, and amplitude summaries.

## Practical Reproduction Notes

- For numerical reproduction, rerun the relevant `build_*` or `run_*` scripts
  with the seed values above.
- For figure verification, use the checked-in `.mat` files and run the listed
  `draw_*` scripts; this avoids long reruns while preserving the exact data used
  by the panels.
- When changing a seed range, also update the saved metadata and this document
  so that repeated simulations remain auditable.
