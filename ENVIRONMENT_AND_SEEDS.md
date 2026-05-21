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
MATLAB Test R2024b
Optimization Toolbox R2024b
Parallel Computing Toolbox R2024b
Signal Processing Toolbox R2024b
Statistics and Machine Learning Toolbox R2024b
Symbolic Math Toolbox R2024b
```

The core `FMAM_ODE` solver only requires MATLAB. Specific workflows use the
following additional components:

- Symbolic Math Toolbox: `build_symbolic_derivatives.m` and scripts that
  regenerate analytic derivative caches.
- Optimization Toolbox: `fmincon` baseline refinement.
- Parallel Computing Toolbox: `parfor` workflows, including repeated SDE
  simulations and large baseline/parameter sweeps.
- Statistics and Machine Learning Toolbox: sampling, distribution summaries, and
  statistical helper functions used by selected plotting or analysis scripts.
- Signal Processing Toolbox: PSD-related analysis in stochastic circadian
  workflows.
- MATLAB Test: repository regression tests under `test/`.
- MatCont 7.6: included in `MatCont7p6/` for workflows using the MATCONT
  periodic-orbit backend.

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
| Baseline comparison | `RLT_circuit/parameter_inference/baseline_comparison/run_baseline_comparison.m` uses `randomSeeds = 1:100`; individual methods call `rng(seed, 'twister')`. |
| Baseline comparison configuration helper | `RLT_circuit/parameter_inference/baseline_comparison/build_config.m` records a default `config.randomSeeds = 0:9`; workflow scripts may override it. |
| Network generation | `network_modulatability/+networkexp/build_network_matrix.m` temporarily sets `rng(seed, 'twister')` and restores the previous RNG state when finished. |
| Network single/ergodic experiments | `network_modulatability/+networkexp/derive_seed.m` and experiment settings determine per-case network and perturbation seeds; saved result files include seed metadata. |
| Circuit sensitivity to initial data | `RLT_circuit/parameter_inference/sensitivity_to_init_data/build_init_data_files.m` saves `seed` and `scale` in each generated initial-data file. |
| Circuit sensitivity to measurement noise | `RLT_circuit/parameter_inference/sensitivity_to_measurement_noise/build_noisy_init_data_files.m` saves `seed` and `noiseLevel` in each generated noisy-data file. |
| Plot jitter only | Scripts such as `draw_baseline_comparison.m` and `draw_sensitivity_to_init_data_violin.m` use fixed `rng(1, 'twister')` before visual scatter jitter; this affects marker placement but not numerical results. |

## Practical Reproduction Notes

- For numerical reproduction, rerun the relevant `build_*` or `run_*` scripts
  with the seed values above.
- For figure verification, use the checked-in `.mat` files and run the listed
  `draw_*` scripts; this avoids long reruns while preserving the exact data used
  by the panels.
- When changing a seed range, also update the saved metadata and this document
  so that repeated simulations remain auditable.
