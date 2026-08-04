# Flexible modulators workflow

FMAM modulation for Fig. 3 (base and temperature models). Tutorial example: `modulation_paths.m`. Condition-number diagnostics: [NearBifurcation/WORKFLOW.md](NearBifurcation/WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd flexible_modulators`
- Symbolic Math Toolbox; MatCont optional for `build_fig3b_data.m`
- Archived `data/fig3/fig3_bifurcation_line.mat` (no builder; used by `draw_fig3c/d` and NearBifurcation)

## Run order

**Plot only:**

```matlab
draw_fig3c; draw_fig3d; draw_fig3f
```

**Full rebuild Fig. 3:**

1. `build_init_data`
2. `build_fig3b_data`
3. Fig. 3c: `build_fig3c_iso_period_curve_data` → `build_fig3c_marker_data` (run 3× with `targetAmplitude` = 2.0, 2.5, 3.0) → `draw_fig3c`
4. Fig. 3d: `build_fig3d_iso_amplitude_curve_data` → `build_fig3d_marker_data` (run 3× with `targetPeriod` = 60, 80, 100) → `draw_fig3d`
5. Fig. 3f: `build_fig3f_data` → `draw_fig3f`

**FMAM tutorial only:** `modulation_paths`

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `modulation_paths.m` | FMAM tutorial (temp model) | `data/modulation_paths/` |
| `build_init_data.m` | Initial orbit | `data/common/initData.mat` |
| `build_fig3b_data.m` | Base-model scan | `data/fig3b/fig3b_base_model_data.mat` |
| `build_fig3c_iso_period_curve_data.m` | Iso-period curves | `data/fig3c/curves/*.mat` |
| `build_fig3c_marker_data.m` | Fig. 3c markers | `data/fig3c/markers/*.mat` |
| `draw_fig3c.m` | Iso-period panel | `fig3c.png` → Fig. 3c |
| `build_fig3d_iso_amplitude_curve_data.m` | Iso-amplitude curves | `data/fig3d/curves/*.mat` |
| `build_fig3d_marker_data.m` | Fig. 3d markers | `data/fig3d/markers/*.mat` |
| `draw_fig3d.m` | Iso-amplitude panel | `fig3d.png` → Fig. 3d |
| `build_fig3f_data.m` | Direct vs orthogonal paths | `data/fig3f/fig3f_data.mat` |
| `draw_fig3f.m` | Path comparison | `fig3f.png` → Fig. 3f |

No `draw_fig3b.m` in this directory.

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. 3c | `draw_fig3c` |
| Fig. 3d | `draw_fig3d` |
| Fig. 3f | `draw_fig3f` |
| Supplementary Fig. S18 (flexible-modulator condition numbers) | [NearBifurcation/WORKFLOW.md](NearBifurcation/WORKFLOW.md) |
