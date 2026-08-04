# RLT circuit tri-target modulation workflow

FMAM period modulation with two amplitudes held fixed (three modulation targets); ODE vs hardware time-series panels (Extended Data Fig. 5, tri-target panels). Parent: [../WORKFLOW.md](../WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `RLT_circuit/learnedData_ODE.mat` from [../parameter_inference/params_inf.m](../parameter_inference/params_inf.m)
- Symbolic Math Toolbox; `FMAM_ODE.m`, `PO_extract/`
- Archived `period_target_*.mat`, `params_modulation_path.mat`, and `circuit_data/*.txt` for plot-only runs

## Run order

**Plot only (archived `.mat` / `circuit_data/`):**

```matlab
cd RLT_circuit/Modulation_tri_target
draw_params
draw_TS   % edit index_PV, period_multiplier, unit at top of script
```

**Full rebuild:**

1. Ensure `../learnedData_ODE.mat` exists (run `parameter_inference/params_inf` if missing).
2. Edit `multiplier` (`1`, `1.5`, or `2`) in `Orthogonal_period_modulation.m`; run three times → `period_target_1p0x.mat`, `period_target_1p5x.mat`, `period_target_2p0x.mat`.
3. Set `needPath = true` and run once → `params_modulation_path.mat` (required by `draw_params`; default is `false`).
4. `draw_params`, `draw_TS` (set `period_multiplier` and `index_PV` per panel).

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `Orthogonal_period_modulation.m` | FMAM period×{1,1.5,2} with `varAmp(1:2)` fixed; controls `R_C`, `inv_C_1`, `inv_C_2` | `period_target_*p*x.mat`; optional `params_modulation_path.mat` |
| `draw_params.m` | 3D \((C_1,C_2,R_C)\) and 2D \((C_3,R_C)\) continuation curves + target markers | interactive figures → Extended Data Fig. 5b |
| `draw_TS.m` | ODE line vs circuit scatter for one state / period multiplier | interactive figure → Extended Data Fig. 5c(i),c(iii) |

Controlled parameters: indices `[1, 4, 5]` (`R_C`, `inv_C_1`, `inv_C_2`). Primary variable: state index 1.

## Figure reproduction

| Figure | Scripts |
|---|---|
| Extended Data Fig. 5b (tri-target parameter curves) | `draw_params` |
| Extended Data Fig. 5c(i),c(iii) (baseline and tri-target time series) | `draw_TS` |

`draw_TS` reads measured traces from `circuit_data/` (`1x_10_20.txt`, `1.5x_10_20.txt`, `2x_10_20.txt`). Bi-target counterpart: [../Modulation_bi_target/WORKFLOW.md](../Modulation_bi_target/WORKFLOW.md).
