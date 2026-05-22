# Network modulatability workflow

Edge-perturbation experiments on GRN and FHN network oscillators (Fig. 1, Extended Data Fig. 1, Figs. S1–S10). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd network_modulatability`
- `PO_extract/`; Parallel Computing Toolbox for ergodic scripts (`useParallel=true`)
- Archived `.mat` in `GRN_SW/`, `FHN_ER/`, `Ergodic data/` enable plot-only runs

## Run order

**Plot only:**

```matlab
% Fig. 1 (GRN): set draw.m top — dynamicName='GRN', netName='SW', N=100, n_per=50
draw; vizNet   % optional topology

% Ext. Data Fig. 1 (FHN): same scripts, dynamicName='FHN'
draw_ergodic
draw_permutation_rank_heatmap   % optional
```

PDFs → `temp_fig/`.

**Full rebuild:**

1. Single experiment: `GRN_network` or `FHN_network` → `draw`, `vizNet`
2. Ergodic sweeps: `GRN_network_ergodic`, `FHN_network_ergodic` → `draw_ergodic`

Top-level scripts define `settings` + `modelSpec`; orchestration in `+networkexp/`.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `GRN_network.m` | Single GRN experiment | `GRN_{net}/TS_init_*.mat`, `TS_per_*.mat` |
| `FHN_network.m` | Single FHN experiment | `FHN_{net}/TS_init_*.mat`, `TS_per_*.mat` |
| `GRN_network_ergodic.m` | GRN robustness sweep | `Ergodic data/N = 100/*.mat` |
| `FHN_network_ergodic.m` | FHN robustness sweep | `Ergodic data/N = 200/*.mat` |
| `draw.m` | Phase/bar/histogram panels | Fig. 1 / Ext. Data Fig. 1 |
| `vizNet.m` | Network + perturbation edges | Fig. 1 / Ext. Data Fig. 1 |
| `draw_ergodic.m` | Success-rate heatmaps | Figs. S1–S10 |
| `draw_permutation_rank_heatmap.m` | Source-node rank heatmap | Figs. S1–S10 supplement |
| `+networkexp/run_single_experiment.m` | Single pipeline | (called by `*_network.m`) |
| `+networkexp/run_ergodic_experiment.m` | Monte Carlo sweep | (called by `*_network_ergodic.m`) |

Match `draw_ergodic.m` settings (`dynamicName`, `netName`, `N`, `sourceSequenceIndex`, `weight_per`) to filenames under `Ergodic data/`.

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 1 (GRN) | `GRN_network` → `draw`, `vizNet` |
| Extended Data Fig. 1 (FHN) | `FHN_network` → `draw`, `vizNet` |
| Figs. S1–S10 | `*_network_ergodic` → `draw_ergodic`, `draw_permutation_rank_heatmap` |

Seeds: `settings.randomSeed` (default 1) via `networkexp.derive_seed`. See [ENVIRONMENT_AND_SEEDS.md](../ENVIRONMENT_AND_SEEDS.md).
