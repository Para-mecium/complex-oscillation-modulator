# RLT circuit baseline comparison workflow

Black-box optimizers vs the proposed FMAM continuation method for electronic-circuit parameter inference (Fig. 4 baseline panels). Parent: [../../WORKFLOW.md](../../WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `../initData_circuit.mat` (circuit measurement target) and `../initData_ODE.mat` (FMAM initial guess)
- Parallel Computing Toolbox recommended for `run_baseline_comparison` (`numWorkers = 8`)
- **Plot coupling:** `draw_baseline_comparison` overlays proposed-method runtime/loss scatter from [../sensitivity_to_init_data/](../sensitivity_to_init_data/) at **scale = 1** (`results/scale_1/`). Run that sub-workflow first (at least `build_sensitivity_to_init_data` for scale 1 + `analyze_sensitivity_to_init_data`).

## Run order

**Plot only (archived `results/`):**

```matlab
cd RLT_circuit/parameter_inference/baseline_comparison
draw_baseline_comparison
```

**Full rebuild:**

1. `params_inf` → `results/proposed_method/relative_l2_orbit/proposed_method_result.mat` (single-run proposed baseline).
2. `run_baseline_comparison` → per-method/per-seed `.mat` under `results/<method>/relative_l2_orbit/` and `results/baseline_comparison_summary_relative_l2_orbit.mat` (500×1000 eval budget × 100 seeds × 5 methods; long runtime).
3. Complete [sensitivity_to_init_data](../sensitivity_to_init_data/WORKFLOW.md) through `analyze_sensitivity_to_init_data` for **scale = 1**.
4. `draw_baseline_comparison`.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `params_inf.m` | Single proposed FMAM inference (continuation from `initData_ODE` → circuit targets) | `results/proposed_method/<loss>/proposed_method_result.mat` |
| `run_baseline_comparison.m` | Parallel sweep: uniform, LHS, DE, Sobol, PSO (+ Powell refinement where configured) | `results/<method>/<loss>/baseline_*_seed_*.mat`, `results/baseline_comparison_summary_<loss>.mat` |
| `draw_baseline_comparison.m` | Best-so-far curves, final-loss violins, runtime scatter, orbit/residual TS | interactive figures → Fig. 4 (baseline panels) |

Shared helpers (`build_config.m`, `baseline_*.m`, `run_powell.m`, `evaluate_candidate.m`, …) are called by the runners above; no separate user entry point.

## Cross-folder dependencies

| Consumer | Required upstream data |
|---|---|
| `draw_baseline_comparison` | `results/baseline_comparison_summary_*.mat`; optional `results/proposed_method/*/proposed_method_result.mat` |
| `draw_baseline_comparison` (proposed scatter) | `../sensitivity_to_init_data/results/scale_1/successful_params_stats.mat` and `params_inf_sensitivity_summary.mat` |

Target orbit: `initData_circuit.mat`. Search bounds and active parameters (`R_C`, `inv_C_1`, `inv_C_2`, `inv_C_3`) are set in `build_config.m` / `run_baseline_comparison.m`.

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 4 (baseline comparison) | `draw_baseline_comparison` (after `run_baseline_comparison` + sensitivity scale 1) |

**Default loss:** `relative_l2_orbit` in `params_inf.m`, `run_baseline_comparison.m`, and `draw_baseline_comparison.m` (matches `build_config.m` and sensitivity analyze). Proposed-method scatter losses use the same `config.lossOptions` when overlaying sensitivity scale 1.
