# RLT workflow

Repressilator FMAM period and phase-difference modulation (Fig. S15). Global map: [FIGURE_WORKFLOW_MAP.md](../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd RLT`
- Symbolic Math Toolbox; `FMAM_ODE.m`, `PO_extract/`
- `.mat` outputs are gitignored; run `run_repressilator_fmam` if `figS15_repressilator_data.mat` is missing

## Run order

**Plot only:**

```matlab
draw_figS15
```

**Full rebuild:**

1. `run_repressilator_fmam` → `initData_repressilator.mat`, `period_target_50.mat`, `phase_target_30.mat`, `phase_target_35.mat`, `figS15_repressilator_data.mat`
2. `draw_figS15` → `Fig_S15.png`, `Fig_S15.pdf`

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `System.m` | Repressilator ODE | — |
| `repressilator_forward_orbit.m` | Periodic orbit extraction | — |
| `run_repressilator_fmam.m` | FMAM period + phase modulation | `figS15_repressilator_data.mat` |
| `draw_figS15.m` | Phase-difference panels | Fig. S15 |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. S15 | `run_repressilator_fmam` (data) → `draw_figS15` |

Edit `periodTarget`, `phaseTargets`, or `Parameters0` in `run_repressilator_fmam.m` to change targets.
