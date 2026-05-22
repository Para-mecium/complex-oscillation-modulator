# RLT circuit sensitivity to measurement noise workflow

FMAM inference from synthetically noisy circuit time series at three noise levels (Fig. 4 measurement-noise panel). Parent: [../WORKFLOW.md](../WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `../initData_circuit.mat` (noiseless target orbit) and `../initData_ODE.mat` (FMAM initial state / parameters)
- Signal Processing Toolbox (`plomb` in `build_noisy_init_data_files`)
- Symbolic Math Toolbox (`build_modulation_results`)
- Independent of `sensitivity_to_init_data` results (only shares parent `initData_*.mat`)

## Run order

**Plot only (archived `noisy_init_data_files/`, `modulation_results/`):**

```matlab
cd RLT_circuit/parameter_inference/sensitivity_to_measurement_noise
draw_sensitivity_to_measurement_noise   % edit errorMetric at top: 'parameter' or 'orbit_l2'
```

**Full rebuild:**

1. `build_noisy_init_data_files` → `noisy_init_data_files/noise_level_*/initData_*.mat` (100 seeds × noise levels `[0.01, 0.02, 0.05]`).
2. `build_modulation_results` → `modulation_results/noise_level_*/modulated_to_initData_*.mat` and `modulation_summary.mat` (100 FMAM runs per level; expensive).
3. `draw_sensitivity_to_measurement_noise`.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_noisy_init_data_files.m` | Add Gaussian noise to `initData_circuit`; Lomb period estimate | `noisy_init_data_files/noise_level_*/initData_*.mat` |
| `build_modulation_results.m` | FMAM from ODE baseline to each noisy target | `modulation_results/noise_level_*/modulated_to_*.mat`, `modulation_summary.mat` |
| `draw_sensitivity_to_measurement_noise.m` | Raincloud of parameter distance or orbit L2 vs noise level | `Sensitivity_to_measurement_noise*.png` → Fig. 4 |

`build_modulation_results` initializes FMAM from `initData_ODE.mat` parameters and fits noisy `period` / `varAmp` targets read from each noisy init file.

## Cross-folder dependencies

| Input | Source |
|---|---|
| Noiseless circuit orbit | `../initData_circuit.mat` |
| FMAM initial parameters / Fourier state | `../initData_ODE.mat` |
| Parameter reference (`errorMetric = 'parameter'`) | `../../learnedData_ODE.mat` (from `../params_inf.m`) |
| Orbit reference (`errorMetric = 'orbit_l2'`) | `../initData_circuit.mat` |

No outputs from this folder are consumed by other RLT_circuit sub-workflows.

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 4 (measurement-noise sensitivity) | `draw_sensitivity_to_measurement_noise` |

Edit `noiseLevels` / `selectedNoiseLevels` and `selectedFiles` in the build scripts to change the sweep. Default draw metric is `orbit_l2`.
