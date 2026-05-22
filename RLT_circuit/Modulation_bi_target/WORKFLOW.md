# RLT circuit bi-target modulation workflow

FMAM period modulation with one amplitude held fixed (two modulation targets); ODE vs hardware time-series panels (Extended Data Fig. 3, bi-target panels). Parent: [../WORKFLOW.md](../WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `RLT_circuit/learnedData_ODE.mat` from [../parameter_inference/params_inf.m](../parameter_inference/params_inf.m) (inferred parameters + ODE orbit)
- Symbolic Math Toolbox; `FMAM_ODE.m`, `PO_extract/`
- Archived `period_target_*.mat`, `params_modulation_path.mat`, and `circuit_data/*.txt` for plot-only runs

## Run order

**Plot only (archived `.mat` / `circuit_data/`):**

```matlab
cd RLT_circuit/Modulation_bi_target
draw_params
draw_TS   % edit index_PV, period_multiplier, unit at top of script
```

**Full rebuild:**

1. Ensure `../learnedData_ODE.mat` exists (run `parameter_inference/params_inf` if missing).
2. Edit `multiplier` (`1`, `1.5`, or `2`) in `Orthogonoal_period_modulation.m`; run three times → `period_target_1p0x.mat`, `period_target_1p5x.mat`, `period_target_2p0x.mat`.
3. Set `needPath = true` and run once (any multiplier that reaches the target) → `params_modulation_path.mat` (continuation log for parameter-plane curves).
4. `draw_params`, `draw_TS` (set `period_multiplier` and `index_PV` per panel).

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `Orthogonoal_period_modulation.m` | FMAM period×{1,1.5,2} with `varAmp(1)` fixed; controls `R_C`, `inv_C_3` | `period_target_*p*x.mat`; optional `params_modulation_path.mat` |
| `draw_params.m` | 3D \((C_1,C_2,R_C)\) and 2D \((C_3,R_C)\) continuation curves + target markers | interactive figures → Extended Data Fig. 3 (parameter panels) |
| `draw_TS.m` | ODE line vs circuit scatter for one state / period multiplier | interactive figure → Extended Data Fig. 3 (time series) |

Controlled parameters: indices `[1, 4]` (`R_C`, `inv_C_3`). Primary variable: state index 1.

## Figure reproduction

| Figure | Scripts |
|---|---|
| Extended Data Fig. 3 (bi-target parameter curves) | `draw_params` |
| Extended Data Fig. 3 (bi-target time series) | `draw_TS` |

`draw_TS` reads measured traces from `circuit_data/` (`1x_10.txt`, `1.5x_10.txt`, `2x_10.txt`). Edit `multiplier`, `needPath`, and `M` in `Orthogonoal_period_modulation.m` to change FMAM settings.
