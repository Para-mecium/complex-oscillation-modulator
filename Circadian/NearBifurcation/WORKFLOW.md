# NearBifurcation workflow (Circadian)

Condition-number diagnostics near Hopf bifurcation on the circadian model (Fig. S17). Parent: [../WORKFLOW.md](../WORKFLOW.md). Flexible-modulator variant: [../../flexible_modulators/NearBifurcation/WORKFLOW.md](../../flexible_modulators/NearBifurcation/WORKFLOW.md). Global map: [FIGURE_WORKFLOW_MAP.md](../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `../data/common/initData.mat` from `../build_init_data.m`
- Symbolic Math Toolbox; `+utils/cbrewer2/` for colormap

## Run order

**Plot only:**

```matlab
cd Circadian/NearBifurcation
draw_iso_period_condition
draw_iso_amplitude_condition
draw_iso_maximum_condition
```

**Full rebuild:**

1. `cd Circadian`; `build_init_data`
2. `cd NearBifurcation`; `build_bifurcation_line_data`
3. `build_iso_period_condition_data`; `build_iso_amplitude_condition_data`; `build_iso_maximum_condition_data`
4. `draw_iso_*_condition` (all three)

**Note:** `build_iso_period_condition_data.m` may list only `targetPeriods = [25.0]` while `draw_iso_period_condition.m` expects four files (T = 23.5, 24.0, 24.5, 25.0). Restore the commented four-period config before full iso-period rebuild; archived `data/iso_period/curves/` already has four files.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_bifurcation_line_data.m` | Analytic Hopf curve | `data/bifurcation/circadian_bifurcation_line.mat` |
| `build_iso_period_condition_data.m` | Iso-period + condition logs | `data/iso_period/curves/*.mat` |
| `build_iso_amplitude_condition_data.m` | Iso-amplitude + condition logs | `data/iso_amplitude/curves/*.mat` |
| `build_iso_maximum_condition_data.m` | Iso-maximum + condition logs | `data/iso_maximum/curves/*.mat` |
| `draw_iso_period_condition.m` | Hopf line + colored curves | `iso_period_condition.png` |
| `draw_iso_amplitude_condition.m` | Same layout | `iso_amplitude_condition.png` |
| `draw_iso_maximum_condition.m` | Same layout | `iso_maximum_condition.png` |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. S17 (iso-period) | `draw_iso_period_condition` |
| Fig. S17 (iso-amplitude) | `draw_iso_amplitude_condition` |
| Fig. S17 (iso-maximum) | `draw_iso_maximum_condition` |
