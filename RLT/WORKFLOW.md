# RLT workflow

Repressilator FMAM period and phase-difference modulation (Supplementary Fig. S17). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd RLT`
- Symbolic Math Toolbox; `FMAM_ODE.m`, `PO_extract/`
- `.mat` outputs are gitignored; run `run_repressilator_fmam` if `figS17_repressilator_data.mat` is missing

## Run order

**Plot only:**

```matlab
draw_figS17
```

**Full rebuild:**

1. `run_repressilator_fmam` → `initData_repressilator.mat`, `period_target_50.mat`, `phase_target_30.mat`, `phase_target_35.mat`, `figS17_repressilator_data.mat`
2. `draw_figS17` → `Fig_S17.png`, `Fig_S17.pdf`

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `System.m` | Repressilator ODE | — |
| `repressilator_forward_orbit.m` | Periodic orbit extraction | — |
| `run_repressilator_fmam.m` | FMAM period + phase modulation | `figS17_repressilator_data.mat` |
| `draw_figS17.m` | Phase-difference panels | Supplementary Fig. S17 |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Supplementary Fig. S17 | `run_repressilator_fmam` (data) → `draw_figS17` |

Edit `periodTarget`, `phaseTargets`, or `Parameters0` in `run_repressilator_fmam.m` to change targets.
