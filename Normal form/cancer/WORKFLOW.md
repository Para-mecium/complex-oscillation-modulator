# Cancer normal-form workflow

Closed-form AM/FM normal-form curves and FMAM comparison for the cancer gene-regulatory oscillator (Fig. S14, cancer panels). Shared utilities: `Normal form/+normalform/`. FHN counterpart: `Normal form/FHN/`. Global map: [FIGURE_WORKFLOW_MAP.md](../../FIGURE_WORKFLOW_MAP.md).

## Prerequisites

- Repo root: `addpath(genpath(pwd))`; `cd 'Normal form/cancer'` (or call `ensure_paths` from any function script here)
- Symbolic Math Toolbox (`build_cancer_normal_form_data.m` → `normalform.build_symbolic_two_state_model`)
- `FMAM_ODE.m`, `PO_extract/` for `run_cancer_fmam.m`

## Run order

**Normal-form curves only:**

```matlab
reproduce_cancer_normal_form
% or: data = build_cancer_normal_form_data; plot_cancer_normal_form(data);
```

**FMAM vs normal-form comparison (needs `cancer_fmam_results.mat`):**

```matlab
plot_cancer_fmam_period_invariant      % AM sweep
plot_cancer_fmam_amplitude_invariant   % FM sweep
```

**Full rebuild (Fig. S14 cancer panels):**

1. `reproduce_cancer_normal_form` → normal-form AM/FM curves (interactive figure)
2. `run_cancer_fmam` → `cancer_fmam_results.mat`
3. `plot_cancer_fmam_period_invariant`, `plot_cancer_fmam_amplitude_invariant` → comparison figures

Step 2 is expensive (41-point AM + 31-point FM continuations); run once and reuse the saved `.mat`.

## Scripts

| Script | Purpose | Output / figure |
|---|---|---|
| `System.m` | Cancer ODE with control parameters (`f11`, `f12`, `f21`) | — |
| `private/ensure_paths.m` | Add repo, `+normalform`, and cancer paths | — |
| `build_cancer_normal_form_data.m` | Closed-form AM/FM curves + regression checks | in-memory `data` struct |
| `plot_cancer_normal_form.m` | Normal-form curve panels | interactive figure |
| `reproduce_cancer_normal_form.m` | Build + optional plot (one call) | `data`, optional figure |
| `run_cancer_fmam.m` | FMAM period- and amplitude-invariant sweeps | `cancer_fmam_results.mat` |
| `plot_cancer_fmam_period_invariant.m` | FMAM AM vs normal-form AM | Fig. S14 (period panel) |
| `plot_cancer_fmam_amplitude_invariant.m` | FMAM FM vs normal-form FM | Fig. S14 (amplitude panel) |
| `cancer.m` | Legacy one-off FMAM run (loads `cancer.mat`, old API) | — |

## Figure reproduction

| Figure | Scripts |
|---|---|
| Fig. S14 (cancer normal-form curves) | `reproduce_cancer_normal_form` |
| Fig. S14 (cancer FMAM comparison) | `run_cancer_fmam` → `plot_cancer_fmam_period_invariant`, `plot_cancer_fmam_amplitude_invariant` |

Edit `amScaleList`, `periodScaleList`, or `baseParameters` in `run_cancer_fmam.m` to change FMAM sweep ranges. Normal-form grids are set in `build_cancer_normal_form_data.m` defaults (`am_f11_range`, `fm_f11_range`). Plot scripts pass `epsilon` from saved FMAM baseline params when rebuilding normal-form data.
