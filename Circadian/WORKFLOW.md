# Circadian workflow

FMAM modulation, SDE stochastic response, and figure scripts for the circadian model (Fig. 5, Fig. S13). Condition-number diagnostics: [NearBifurcation/WORKFLOW.md](NearBifurcation/WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- Symbolic Math Toolbox (`build_*`); Parallel Computing + Signal Processing Toolboxes for `build_fig5d_sde_repeat_data.m`
- Seeds: [ENVIRONMENT_AND_SEEDS.md](../ENVIRONMENT_AND_SEEDS.md)

## Run order

**Plot only (archived `data/`):**

```matlab
cd Circadian
draw_fig5b; draw_fig5c; draw_fig5d
draw_fig5d_sde_representative; draw_fig5d_sde_stats
draw_figS2a; draw_figS2b   % S2b needs data/fig5d/curves/
```

**Full rebuild:**

1. `build_init_data`
2. Fig. 5b–c: `build_fig5b_iso_period_curve_data` → `build_fig5b_marker_data` → `draw_fig5b`, `draw_fig5c`
3. Fig. 5d: `build_fig5d_iso_maximum_curve_data` → `build_fig5d_marker_data` → `draw_fig5d`
4. Fig. 5d SDE: `build_fig5d_sde_representative_data` → `build_fig5d_sde_repeat_data` → `draw_fig5d_sde_*`
5. Fig. S13: `build_figS2a_iso_amplitude_curve_data` → `build_figS2a_marker_data` → `draw_figS2a`; `draw_figS2b` (shares fig5d curves)

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_init_data.m` | Initial periodic orbit | `data/common/initData.mat` |
| `build_fig5b_iso_period_curve_data.m` | Iso-period curves | `data/fig5b/curves/*.mat` |
| `build_fig5b_marker_data.m` | Fig. 5b markers | `data/fig5b/markers/*.mat` |
| `draw_fig5b.m` | Parameter plane | `fig5b.png` → Fig. 5b |
| `draw_fig5c.m` | Marker time series | `fig5c.png` → Fig. 5c |
| `build_fig5d_iso_maximum_curve_data.m` | Iso-maximum curves | `data/fig5d/curves/*.mat` |
| `build_fig5d_marker_data.m` | Fig. 5d markers | `data/fig5d/markers/*.mat` |
| `draw_fig5d.m` | Parameter plane | `fig5d.png` → Fig. 5d |
| `build_fig5d_sde_representative_data.m` | Representative SDE | `data/fig5d/sde/fig5d_sde_representative.mat` |
| `build_fig5d_sde_repeat_data.m` | 100× SDE ensemble | `data/fig5d/sde/fig5d_sde_repeat.mat` |
| `draw_fig5d_sde_representative.m` | SDE trajectories | Fig. 5d SDE panel |
| `draw_fig5d_sde_stats.m` | PSD + distribution | Fig. 5d SDE stats |
| `build_figS2a_iso_amplitude_curve_data.m` | Iso-amplitude curves | `data/figS2a/curves/*.mat` |
| `build_figS2a_marker_data.m` | Fig. S13a markers | `data/figS2a/markers/*.mat` |
| `draw_figS2a.m` | Iso-amplitude panel | `figS2a.png` → Fig. S13 |
| `draw_figS2b.m` | Iso-maximum panel | `figS2b.png` → Fig. S13 |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 5b–c | `draw_fig5b`, `draw_fig5c` |
| Fig. 5d (deterministic) | `draw_fig5d` |
| Fig. 5d (SDE) | `draw_fig5d_sde_representative`, `draw_fig5d_sde_stats` |
| Fig. S13 | `draw_figS2a`, `draw_figS2b` |
| Fig. S17 | [NearBifurcation/WORKFLOW.md](NearBifurcation/WORKFLOW.md) |
