# Longevity workflow

FMAM parameter modulation and SDE validation for the yeast longevity oscillator (Extended Data Fig. 3 and Supplementary Fig. S13). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd Longevity`
- Symbolic Math Toolbox for `build_*_data.m`; Statistics Toolbox for `draw_ed_fig3d.m`
- Archived `.mat` files in `Longevity/` enable plot-only runs

## Run order

**Plot only:**

```matlab
draw_ed_fig3b; draw_ed_fig3c; draw_ed_fig3d; draw_ed_fig3e
draw_figS13a
```

**Full rebuild:**

1. `build_init_data`
2. `build_ed_fig3b_data` → `build_ed_fig3c_data` → `build_ed_fig3d_data` → `draw_ed_fig3b`–`draw_ed_fig3e`
3. Supplementary Fig. S13: `build_figS13a_data` → `draw_figS13a` (shares `initData.mat`)

If SDE `.mat` files are missing, `draw_ed_fig3c` and `draw_ed_fig3d` run their matching builders automatically.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_init_data.m` | Initial orbit | `initData.mat` |
| `build_ed_fig3b_data.m` | α modulation path | `alpha_modulation_path.mat`, `alpha_target_data.mat` |
| `build_ed_fig3c_data.m` | Representative SDE | `sde_representative.mat` |
| `build_ed_fig3d_data.m` | 100× ensemble SDE | `sde_distribution.mat` |
| `draw_ed_fig3b.m` | Modulation panels | Extended Data Fig. 3b |
| `draw_ed_fig3c.m` | SDE time series | Extended Data Fig. 3c |
| `draw_ed_fig3d.m` | Min S/H histograms | Extended Data Fig. 3d |
| `draw_ed_fig3e.m` | Ensemble bands | Extended Data Fig. 3c (right-column summary) |
| `build_figS13a_data.m` | β modulation path | `beta_modulation_path.mat`, `beta_target_data.mat` |
| `draw_figS13a.m` | β modulation panels | Supplementary Fig. S13 |
| `phase_plot.m` | Legacy interactive phase plot | Extended Data Fig. 3b helper (needs workspace FMAM output) |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Extended Data Fig. 3b | `draw_ed_fig3b` |
| Extended Data Fig. 3c–d | `draw_ed_fig3c`, `draw_ed_fig3d`, `draw_ed_fig3e` (`draw_ed_fig3e` supplies the right-column summary in c) |
| Supplementary Fig. S13 | `draw_figS13a` |

SDE seeds: `build_ed_fig3c_data` uses `seed=1`; `build_ed_fig3d_data` uses `1001:1100`. See [ENVIRONMENT_AND_SEEDS.md](../ENVIRONMENT_AND_SEEDS.md).
