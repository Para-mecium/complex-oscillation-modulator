# NearBifurcation workflow (flexible modulators)

Condition-number diagnostics along iso-period/iso-amplitude FMAM curves on the base flexible modulator (Supplementary Fig. S18). Parent: [../WORKFLOW.md](../WORKFLOW.md). Iso-maximum variant is in `Circadian/NearBifurcation/`. Global map: [FIGURE_WORKFLOW_MAP.md](../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`
- `../data/common/initData.mat` from `../build_init_data.m`
- Archived `../data/fig3/fig3_bifurcation_line.mat` (Hopf/LP lines; no builder here)
- Symbolic Math Toolbox; `+utils/cbrewer2/` for colormap

Builders use `needLog=true`, `conditionStopRcond=1e-12` (stricter than Fig. 3c/3d builders). Data is independent of `data/fig3c/` and `data/fig3d/`.

## Run order

**Plot only:**

```matlab
cd flexible_modulators/NearBifurcation
draw_iso_period_condition
draw_iso_amplitude_condition
```

**Full rebuild:**

1. `cd flexible_modulators`; `build_init_data`
2. `cd NearBifurcation`; `build_iso_period_condition_data`; `build_iso_amplitude_condition_data`
3. `draw_iso_period_condition`; `draw_iso_amplitude_condition`

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `build_iso_period_condition_data.m` | Iso-period + condition logs | `data/iso_period/curves/iso_period_T*_condition.mat` |
| `build_iso_amplitude_condition_data.m` | Iso-amplitude + condition logs | `data/iso_amplitude/curves/iso_amplitude_A*_condition.mat` |
| `draw_iso_period_condition.m` | log10(condition) coloring | `iso_period_condition.png` → Supplementary Fig. S18b |
| `draw_iso_amplitude_condition.m` | log10(condition) coloring | `iso_amplitude_condition.png` → Supplementary Fig. S18a |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Supplementary Fig. S18b (iso-period) | `draw_iso_period_condition` |
| Supplementary Fig. S18a (iso-amplitude) | `draw_iso_amplitude_condition` |

Targets: periods `T ∈ {50,60,70,80,90}`; amplitudes `A = 1.2:0.3:3.6`.
