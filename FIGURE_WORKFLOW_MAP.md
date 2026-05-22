# Figure Workflow Map

This file maps the manuscript figure groups to the repository directories and
scripts. The map distinguishes data generation, inference or simulation, and
plotting steps where those steps are separate.

Run scripts from the repository root after:

```matlab
addpath(genpath(pwd))
```

## Main Figure Groups

| Figure group | Repository contents | Data generation / inference | Plotting scripts |
|---|---|---|---|
| Fig. 1, GRN network modulatability | `network_modulatability/` | `GRN_network.m`, `+networkexp/run_single_experiment.m`, `+networkexp/run_ergodic_experiment.m` | `draw.m`, `draw_ergodic.m`, `vizNet.m` |
| Fig. 2, longevity-oriented modulation | `Longevity/` | `build_init_data.m`, `build_etd_fig2b_data.m`, `build_etd_fig2c_data.m`, `build_etd_fig2d_data.m` | `draw_etd_fig2b.m`, `draw_etd_fig2c.m`, `draw_etd_fig2d.m`, `draw_etd_fig2e.m`, `phase_plot.m` |
| Fig. 3, flexible modulator | `flexible_modulators/` | `build_init_data.m`, `build_fig3b_data.m`, `build_fig3c_*_data.m`, `build_fig3d_*_data.m`, `build_fig3f_data.m` | `draw_fig3c.m`, `draw_fig3d.m`, `draw_fig3f.m` |
| Fig. 4, electronic circuit parameter inference | `RLT_circuit/parameter_inference/` | `params_inf.m`, `baseline_comparison/params_inf.m`, `baseline_comparison/run_baseline_comparison.m`, `sensitivity_to_init_data/build_*.m`, `sensitivity_to_measurement_noise/build_*.m` | `draw_TS.m`, `draw_params.m`, `baseline_comparison/draw_baseline_comparison.m`, `sensitivity_to_*/draw_*.m` |
| Fig. 5, circadian modulation and stochastic response | `Circadian/` | `build_init_data.m`, `build_fig5b_*_data.m`, `build_fig5d_*_data.m`, `build_fig5d_sde_representative_data.m`, `build_fig5d_sde_repeat_data.m` | `draw_fig5b.m`, `draw_fig5c.m`, `draw_fig5d.m`, `draw_fig5d_sde_representative.m`, `draw_fig5d_sde_stats.m` |

## Extended Data Figure Groups

| Figure group | Repository contents | Data generation / inference | Plotting scripts |
|---|---|---|---|
| Extended Data Fig. 1, FHN network modulatability | `network_modulatability/` | `FHN_network.m`, `+networkexp/run_single_experiment.m`, `+networkexp/run_ergodic_experiment.m` | `draw.m`, `draw_ergodic.m`, `vizNet.m` |
| Extended Data Fig. 2, engineered longevity oscillator | `Longevity/` | `build_init_data.m`, `build_etd_fig2b_data.m`, `build_etd_fig2c_data.m`, `build_etd_fig2d_data.m` | `draw_etd_fig2b.m`, `draw_etd_fig2c.m`, `draw_etd_fig2d.m`, `draw_etd_fig2e.m`, `phase_plot.m` |
| Extended Data Fig. 3, RLT circuit orthogonal modulation | `RLT_circuit/Modulation_bi_target/`, `RLT_circuit/Modulation_tri_target/` | `Orthogonoal_period_modulation.m`, `Orthogonal_period_modulation.m` | `draw_params.m`, `draw_TS.m` |
| Extended Data Fig. 4, complexity sweep | `complexity_analyses/` | `run_fhn_complexity_analysis.m` | `plot_fhn_complexity_results.m`, `reproduce_fhn_complexity_sweeps.m` |
| Extended Data Fig. 5, inferred-parameter distributions | `RLT_circuit/parameter_inference/sensitivity_to_init_data/` | `build_init_data_files.m`, `build_sensitivity_to_init_data.m`, `analyze_sensitivity_to_init_data.m` | `draw_sensitivity_to_init_data.m` |

## Principal Supplementary Figure Groups

| Figure group | Repository contents | Data generation / inference | Plotting scripts |
|---|---|---|---|
| Figs. S1-S10, network robustness sweeps | `network_modulatability/` | `GRN_network.m`, `FHN_network.m`, `GRN_network_ergodic.m`, `FHN_network_ergodic.m`, `+networkexp/*.m` | `draw.m`, `draw_ergodic.m`, `draw_permutation_rank_heatmap.m` |
| Fig. S11, longevity supplementary modulation | `Longevity/` | `build_init_data.m`, `build_figS1a_data.m` | `draw_figS1a.m` |
| Fig. S13, circadian supplementary modulation | `Circadian/` | `build_figS2a_*_data.m`, `build_fig5d_*_data.m` | `draw_figS2a.m`, `draw_figS2b.m` |
| Fig. S14, normal-form comparison | `Normal form/FHN/`, `Normal form/cancer/` | `run_fhn_fmam.m`, `run_cancer_fmam.m`, `build_fhn_normal_form_data.m`, `build_cancer_normal_form_data.m` | `plot_fhn_normal_form.m`, `plot_cancer_normal_form.m`, `plot_*_fmam_amplitude_invariant.m`, `plot_*_fmam_period_invariant.m`, `reproduce_fhn_normal_form.m`, `reproduce_cancer_normal_form.m` |
| Fig. S15, repressilator phase-difference example | `RLT/` | `run_repressilator_fmam.m` | `draw_figS15.m` |
| Figs. S16-S17, condition-number diagnostics | `flexible_modulators/NearBifurcation/`, `Circadian/NearBifurcation/` | `build_iso_period_condition_data.m`, `build_iso_amplitude_condition_data.m`, `build_iso_maximum_condition_data.m`, `build_bifurcation_line_data.m` | `draw_iso_period_condition.m`, `draw_iso_amplitude_condition.m`, `draw_iso_maximum_condition.m` |

## Reading The Map

- If a figure can be verified from archived `.mat` data, run the draw script
  first.
- If regenerated numerical data are required, run the listed builder or
  inference script before the draw script.
- For stochastic or repeated simulations, the random seeds are documented in
  `ENVIRONMENT_AND_SEEDS.md` and usually saved into the corresponding `.mat`
  output.
- For computationally expensive sweeps, archived intermediate data are part of
  the reproducibility package; full regeneration may require substantially more
  runtime than plotting from the archived files.
