# FHN normal-form workflow

Closed-form AM/FM normal-form curves and FMAM comparison for the FitzHugh–Nagumo oscillator (Supplementary Fig. S16, FHN data). Shared utilities: `Normal form/+normalform/`. Cancer counterpart: `Normal form/cancer/`. Global map: [FIGURE_WORKFLOW_MAP.md](../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd 'Normal form/FHN'` (or call `ensure_paths` from any script here)
- Symbolic Math Toolbox (`build_fhn_normal_form_data.m` → `normalform.build_symbolic_two_state_model`)
- `FMAM_ODE.m`, `PO_extract/` for `run_fhn_fmam.m`

## Run order

**Normal-form curves only:**

```matlab
reproduce_fhn_normal_form
% or: data = build_fhn_normal_form_data; plot_fhn_normal_form(data);
```

**FMAM vs normal-form comparison (needs `fhn_fmam_results.mat`):**

```matlab
plot_fhn_fmam_period_invariant      % AM sweep
plot_fhn_fmam_amplitude_invariant   % FM sweep
```

**Full rebuild (Supplementary Fig. S16, FHN data):**

1. `reproduce_fhn_normal_form` → normal-form AM/FM curves (interactive figure)
2. `run_fhn_fmam` → `fhn_fmam_results.mat`
3. `plot_fhn_fmam_period_invariant`, `plot_fhn_fmam_amplitude_invariant` → comparison figures

Step 2 is expensive (101-point AM + 101-point FM continuations); run once and reuse the saved `.mat`.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `System.m` | FHN ODE with control parameters | — |
| `ensure_paths.m` | Add repo, `+normalform`, and FHN paths | — |
| `build_fhn_normal_form_data.m` | Closed-form AM/FM curves + regression checks | in-memory `data` struct |
| `plot_fhn_normal_form.m` | Normal-form curve panels | interactive figure |
| `reproduce_fhn_normal_form.m` | Build + optional plot (one call) | `data`, optional figure |
| `run_fhn_fmam.m` | FMAM period- and amplitude-invariant sweeps | `fhn_fmam_results.mat` |
| `plot_fhn_fmam_period_invariant.m` | FMAM AM vs normal-form AM | Supplementary Fig. S16 (period-invariant data) |
| `plot_fhn_fmam_amplitude_invariant.m` | FMAM FM vs normal-form FM | Supplementary Fig. S16 (amplitude-invariant data) |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Supplementary Fig. S16 (FHN normal-form curves) | `reproduce_fhn_normal_form` |
| Supplementary Fig. S16 (FHN FMAM comparison) | `run_fhn_fmam` → `plot_fhn_fmam_period_invariant`, `plot_fhn_fmam_amplitude_invariant` |

Edit `amScaleList`, `periodScaleList`, or `baseParameters` in `run_fhn_fmam.m` to change FMAM sweep ranges. Normal-form grids are set in `build_fhn_normal_form_data.m` defaults (`am_f11_range`, `fm_f11_range`).
