# Longevity workflow

FMAM parameter modulation and SDE validation for the yeast longevity oscillator (Fig. 2, Extended Data Fig. 2, Fig. S11). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd Longevity`
- Symbolic Math Toolbox for `build_*_data.m`; Statistics Toolbox for `draw_etd_fig2d.m`
- Archived `.mat` files in `Longevity/` enable plot-only runs

## Run order

**Plot only:**

```matlab
draw_etd_fig2b; draw_etd_fig2c; draw_etd_fig2d; draw_etd_fig2e
draw_figS1a
```

**Full rebuild:**

1. `build_init_data`
2. `build_etd_fig2b_data` → `build_etd_fig2c_data` → `build_etd_fig2d_data` → `draw_etd_fig2b`–`draw_etd_fig2e`
3. Fig. S11: `build_figS1a_data` → `draw_figS1a` (shares `initData.mat`)

If SDE `.mat` files are missing, run builders manually (`draw_etd_fig2c/d` fallback names do not match builders).

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_init_data.m` | Initial orbit | `initData.mat` |
| `build_etd_fig2b_data.m` | α modulation path | `alpha_modulation_path.mat`, `alpha_target_data.mat` |
| `build_etd_fig2c_data.m` | Representative SDE | `sde_representative.mat` |
| `build_etd_fig2d_data.m` | 100× ensemble SDE | `sde_distribution.mat` |
| `draw_etd_fig2b.m` | Modulation panels | Fig. 2 / Ext. Data Fig. 2b |
| `draw_etd_fig2c.m` | SDE time series | Ext. Data Fig. 2c |
| `draw_etd_fig2d.m` | Min S/H histograms | Ext. Data Fig. 2d |
| `draw_etd_fig2e.m` | Ensemble bands | Ext. Data Fig. 2e |
| `build_figS1a_data.m` | β modulation path | `beta_modulation_path.mat`, `beta_target_data.mat` |
| `draw_figS1a.m` | β modulation panels | Fig. S11 |
| `phase_plot.m` | Legacy interactive phase plot | Fig. 2 (needs workspace FMAM output) |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 2 / Ext. Data Fig. 2b | `draw_etd_fig2b` |
| Ext. Data Fig. 2c–d | `draw_etd_fig2c`, `draw_etd_fig2d`, `draw_etd_fig2e` (e is merged into the right column of c)|
| Fig. S11 | `draw_figS1a` |

SDE seeds: `build_etd_fig2c_data` uses `seed=1`; `build_etd_fig2d_data` uses `1001:1100`. See [ENVIRONMENT_AND_SEEDS.md](../ENVIRONMENT_AND_SEEDS.md).
