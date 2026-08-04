# RLT circuit sensitivity to initial data workflow

FMAM parameter inference repeated over random initial guesses at five sampling-region scales (Fig. 4 sensitivity panels; Extended Data Fig. 4). Parent: [../WORKFLOW.md](../WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `../initData_circuit.mat` (modulation targets) and `../initData_ODE.mat` (parameter template / `y0`)
- Symbolic Math Toolbox (`build_sensitivity_to_init_data`)
- **Downstream coupling:** [../baseline_comparison/draw_baseline_comparison.m](../baseline_comparison/draw_baseline_comparison.m) reads **scale = 1** results from `results/scale_1/`

## Run order

**Plot only (archived `init_data_files/`, `results/`):**

```matlab
cd RLT_circuit/parameter_inference/sensitivity_to_init_data
draw_sensitivity_to_init_data
draw_sensitivity_to_init_data_violin
draw_sensitivity_to_init_data_success_rate
draw_init_data_outlier_timeseries
```

**Full rebuild:**

1. `build_init_data_files` → `init_data_files/scale_*/initData_*.mat` (100 accepted guesses × scales `[1, 0.5, 0.75, 1.25, 1.5]`).
2. `build_sensitivity_to_init_data` → `results/scale_*/params_inf_sensitivity_summary.mat` (100 FMAM runs per scale; expensive).
3. `analyze_sensitivity_to_init_data` → `results/scale_*/successful_params_stats.mat`.
4. Draw scripts (step 1 plot-only block).

Optional forward-perturbation branch (not in figure map): copy or symlink `results/scale_1/successful_params_stats.mat` to `successful_params_stats.mat` in this folder, then `forward_perturb_params` → `analyze_forward_perturb_params` → `draw_forward_perturb_params`.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_init_data_files.m` | Sample perturbed \((R_C,C_1,C_2,C_3)\) around ODE baseline; forward-orbit filter | `init_data_files/scale_*/initData_*.mat` |
| `build_sensitivity_to_init_data.m` | FMAM continuation from each init guess to circuit targets | `results/scale_*/params_inf_sensitivity_summary.mat` |
| `analyze_sensitivity_to_init_data.m` | Aggregate successes; orbit L2, physical params, correlations | `results/scale_*/successful_params_stats.mat` |
| `draw_sensitivity_to_init_data.m` | Histograms, boxplots, scatter (scale = 1) | `results/scale_1/successful_params_*.png` → Extended Data Fig. 4 |
| `draw_sensitivity_to_init_data_violin.m` | Raincloud over all scales (`errorMetric`: `parameter` or `orbit_l2`) | `results/Sensitivity_to_initial_data_violin.png` → Fig. 4 |
| `draw_sensitivity_to_init_data_success_rate.m` | Continuation success rate vs scale | `results/Sensitivity_to_initial_data_success_rate.png` |
| `draw_init_data_outlier_timeseries.m` | Outlier orbit time series (uses `learnedData_ODE.mat` reference) | `results/outlier_timeseries/*.png` |
| `forward_perturb_params.m` | Forward ODE perturbation from median inferred params | `forward_perturbation_summary.mat` |
| `analyze_forward_perturb_params.m` | Summarize perturbation runs | `forward_perturbation_stats.mat` |
| `draw_forward_perturb_params.m` | Perturbation heatmap / relative-change plots | `forward_perturb_*.png` |

## Cross-folder dependencies

| Output | Used by |
|---|---|
| `results/scale_1/successful_params_stats.mat` | `baseline_comparison/draw_baseline_comparison` (proposed runtime/loss overlay) |
| `results/scale_1/params_inf_sensitivity_summary.mat` | `baseline_comparison/draw_baseline_comparison` (failed-run runtimes) |
| `../learnedData_ODE.mat` | `draw_sensitivity_to_init_data_violin` (`parameter` metric), `draw_init_data_outlier_timeseries` |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 4 (initial-data sensitivity raincloud) | `draw_sensitivity_to_init_data_violin` |
| Fig. 4 (success rate vs scale) | `draw_sensitivity_to_init_data_success_rate` |
| Extended Data Fig. 4 | `draw_sensitivity_to_init_data` |
| Outlier diagnostics | `draw_init_data_outlier_timeseries` |

Scale directories use the pattern `scale_<token>` (e.g. `scale_1`, `scale_0p5`).
