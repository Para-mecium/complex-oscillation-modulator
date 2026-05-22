# Extended Data Fig. 4 — FHN complexity sweep

Systematic FMAM runtime sweeps on coupled FHN networks. Parameter details: [README.md](README.md). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- Depends on `FMAM_ODE.m`, `network_modulatability/+networkexp/`, `PO_extract/`
- Archived data: `data/base_N20_M10/`

## Run order

**Plot from archive (recommended):**

```matlab
reproduce_fhn_complexity_sweeps   % paper layout, all five sweeps
```

Optional single-sweep exploration: edit `User settings` in `plot_fhn_complexity_results.m`, then run it.

**Full rebuild:**

1. Edit `User settings` in `run_fhn_complexity_analysis.m`
2. `run_fhn_complexity_analysis` → writes `data/<baseline_label>/`
3. `reproduce_fhn_complexity_sweeps` (update paths if baseline label changed)

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `run_fhn_complexity_analysis.m` | Batch benchmark sweeps | `data/base_N20_M10/*_sweep/*.mat` |
| `benchmark_fhn_case.m` | Single-case FMAM pipeline | (called by runner) |
| `build_heterogeneous_fhn_network.m` | FHN network model | (called by benchmark) |
| `reproduce_fhn_complexity_sweeps.m` | Paper figures | `figures/base_N20_M10/paper_sweeps/fig_S11–S15_*.png` |
| `plot_fhn_complexity_results.m` | Generic runtime curves | `figures/<baseline>/<sweep>/` |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Extended Data Fig. 4 (system dimension) | `reproduce_fhn_complexity_sweeps` → `fig_S11_system_dimension.png` |
| Extended Data Fig. 4 (truncation order) | → `fig_S12_truncation_order.png` |
| Extended Data Fig. 4 (target number) | → `fig_S13_target_number.png` |
| Extended Data Fig. 4 (Newton tolerance) | → `fig_S14_newton_tolerance.png` |
| Extended Data Fig. 4 (step cap) | → `fig_S15_step_cap.png` |

Each PNG has two subpanels (total runtime; runtime per Newton iteration). Default sweep: 32 cases (~minutes to hours to regenerate; N=100 cases dominate).
