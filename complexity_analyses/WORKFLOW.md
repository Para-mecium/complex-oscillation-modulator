# Extended Data Fig. 2 — FHN complexity sweep

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
| `reproduce_fhn_complexity_sweeps.m` | Paper figures | `figures/base_N20_M10/paper_sweeps/ed_fig2a–2e_*.png` |
| `plot_fhn_complexity_results.m` | Generic runtime curves | `figures/<baseline>/<sweep>/` |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Extended Data Fig. 2 | `run_fhn_complexity_analysis` optionally regenerates the five sweep datasets; `reproduce_fhn_complexity_sweeps` generates five two-panel component PNGs. The current composite figure `fig_extend_complexity.pdf` is assembled separately from these components. |

Each PNG has two subpanels (total runtime; runtime per Newton iteration). Default sweep: 32 cases (~minutes to hours to regenerate; N=100 cases dominate).
