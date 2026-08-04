# Circadian workflow

FMAM modulation, SDE stochastic response, and figure scripts for the circadian model (Fig. 5 and Supplementary Fig. S15). Condition-number diagnostics (Supplementary Fig. S19): [NearBifurcation/WORKFLOW.md](NearBifurcation/WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

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
draw_figS15a; draw_figS15b   % S15b needs data/fig5d/curves/
```

**Full rebuild:**

1. `build_init_data`
2. Fig. 5b–c: `build_fig5b_iso_period_curve_data` → `build_fig5b_marker_data` → `draw_fig5b`, `draw_fig5c`
3. Fig. 5d: `build_fig5d_iso_maximum_curve_data` → `build_fig5d_marker_data` → `draw_fig5d`
4. Fig. 5d SDE: `build_fig5d_sde_representative_data` → `build_fig5d_sde_repeat_data` → `draw_fig5d_sde_*`
5. Supplementary Fig. S15: `build_figS15a_iso_amplitude_curve_data` → `build_figS15a_marker_data` → `draw_figS15a`; `draw_figS15b` (shares fig5d curves)

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
| `build_figS15a_iso_amplitude_curve_data.m` | Iso-amplitude curves | `data/figS15a/curves/*.mat` |
| `build_figS15a_marker_data.m` | Supplementary Fig. S15a markers | `data/figS15a/markers/*.mat` |
| `draw_figS15a.m` | Iso-amplitude panel | `figS15a.png` → Supplementary Fig. S15a |
| `draw_figS15b.m` | Iso-maximum panel | `figS15b.png` → Supplementary Fig. S15b |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 5b–c | `draw_fig5b`, `draw_fig5c` |
| Fig. 5d (deterministic) | `draw_fig5d` |
| Fig. 5d (SDE) | `draw_fig5d_sde_representative`, `draw_fig5d_sde_stats` |
| Supplementary Fig. S15 | `draw_figS15a`, `draw_figS15b` |
| Supplementary Fig. S19 | [NearBifurcation/WORKFLOW.md](NearBifurcation/WORKFLOW.md) |
