# Figure Workflow Map

This file maps the manuscript figure groups to the repository directories,
scripts, archived inputs, and expected outputs. The map distinguishes data
generation, inference or simulation, and plotting steps where those steps are
separate.

Run scripts from the repository root after:

```matlab
addpath(genpath(pwd))
```

## Main Figure Groups

| Figure group | Repository contents | Data generation / inference | Plotting scripts | Archived inputs and expected outputs |
|---|---|---|---|---|
| Fig. 1, GRN network modulatability | `network_modulatability/` | `GRN_network.m`, `+networkexp/run_single_experiment.m`, `+networkexp/run_ergodic_experiment.m` | `draw.m`, `draw_ergodic.m`, `vizNet.m` | Inputs and generated `.mat` files under `GRN_*` and `Ergodic data/`; output panels under `temp_fig/` and `Figures/Fig 1 (GRN network)/` |
| Fig. 2, longevity-oriented modulation | `Longevity/` | `build_init_data.m`, `build_etd_fig2b_data.m`, `build_etd_fig2c_data.m`, `build_etd_fig2d_data.m` | `draw_etd_fig2b.m`, `draw_etd_fig2c.m`, `draw_etd_fig2d.m`, `draw_etd_fig2e.m`, `phase_plot.m` | Inputs include `initData.mat`, `alpha_target_data.mat`, `beta_target_data.mat`, and `sde_*.mat`; expected output is the longevity figure panels |
| Fig. 3, flexible modulator | `flexible_modulators/` | `build_init_data.m`, `build_fig3b_data.m`, `build_fig3c_*_data.m`, `build_fig3d_*_data.m`, `build_fig3f_data.m`, `modulation_paths.m` | `draw_fig3c.m`, `draw_fig3d.m`, `draw_fig3f.m` | Inputs and generated data under `data/common/`, `data/fig3*/`, and `data/modulation_paths/`; expected outputs include `fig3*.png` |
| Fig. 4, electronic circuit parameter inference | `RLT_circuit/parameter_inference/` | `params_inf.m`, `baseline_comparison/params_inf.m`, `baseline_comparison/run_baseline_comparison.m`, `sensitivity_to_init_data/build_*.m`, `sensitivity_to_measurement_noise/build_*.m` | `draw_TS.m`, `draw_params.m`, `baseline_comparison/draw_baseline_comparison.m`, `sensitivity_to_*/draw_*.m` | Inputs include `initData_ODE.mat`, `initData_circuit.mat`, `learnedData_ODE.mat`, baseline `results/`, and sensitivity result folders; expected outputs include baseline and sensitivity panels |
| Fig. 5, circadian modulation and stochastic response | `Circadian/` | `build_init_data.m`, `build_fig5b_*_data.m`, `build_fig5d_*_data.m`, `build_fig5d_sde_representative_data.m`, `build_fig5d_sde_repeat_data.m` | `draw_fig5b.m`, `draw_fig5c.m`, `draw_fig5d.m`, `draw_fig5d_sde_representative.m`, `draw_fig5d_sde_stats.m` | Inputs and generated data under `data/common/`, `data/fig5b/`, `data/fig5d/`; expected outputs include `fig5*.png` and SDE PSD/distribution panels |

## Principal Supplementary Figure Groups

| Figure group | Repository contents | Data generation / inference | Plotting scripts | Archived inputs and expected outputs |
|---|---|---|---|---|
| Extended network/FHN modulatability panels | `network_modulatability/` | `FHN_network.m`, `+networkexp/*.m` | `draw.m`, `draw_ergodic.m`, `draw_permutation_rank_heatmap.m` | Inputs and outputs under `FHN_*`, `Ergodic data/`, and `temp_fig/` |
| Longevity SDE supplementary panels | `Longevity/` | `build_etd_fig2c_data.m`, `build_etd_fig2d_data.m`, `build_figS1a_data.m` | `draw_etd_fig2c.m`, `draw_etd_fig2d.m`, `draw_etd_fig2e.m`, `draw_figS1a.m` | Inputs include `sde_representative.mat` and `sde_distribution.mat`; expected outputs are representative traces, histograms, and percentile/mean bands |
| Circadian supplementary and near-bifurcation panels | `Circadian/`, `Circadian/NearBifurcation/` | `build_figS2a_*_data.m`, `NearBifurcation/build_*_condition_data.m`, `NearBifurcation/build_bifurcation_line_data.m` | `draw_figS2a.m`, `draw_figS2b.m`, `NearBifurcation/draw_*_condition.m` | Inputs under `Circadian/data/figS2a/` and `Circadian/NearBifurcation/data/`; expected outputs include `figS2*.png` and condition-number panels |
| Flexible-modulator near-bifurcation panels | `flexible_modulators/NearBifurcation/` | `build_iso_period_condition_data.m`, `build_iso_amplitude_condition_data.m` | `draw_iso_period_condition.m`, `draw_iso_amplitude_condition.m` | Inputs under `flexible_modulators/NearBifurcation/data/`; expected outputs include `iso_*_condition.png` |
| Complexity analysis panels | `complexity_analyses/` | `run_fhn_complexity_analysis.m` | `plot_fhn_complexity_results.m`, `reproduce_fhn_complexity_sweeps.m` | Inputs under `complexity_analyses/data/`; expected outputs under `complexity_analyses/figures/` |
| Parameter-inference robustness panels | `RLT_circuit/parameter_inference/sensitivity_to_init_data/`, `sensitivity_to_measurement_noise/` | `build_init_data_files.m`, `build_sensitivity_to_init_data.m`, `build_noisy_init_data_files.m`, `build_modulation_results.m` | `draw_sensitivity_to_init_data*.m`, `draw_sensitivity_to_measurement_noise.m`, `draw_measurement_noise_param_histograms.m` | Inputs and summaries under `init_data_files/`, `results/`, `noisy_init_data_files/`, and `modulation_results/`; expected outputs are success-rate, violin, histogram, and outlier panels |
| Baseline comparison panels | `RLT_circuit/parameter_inference/baseline_comparison/` | `run_baseline_comparison.m`, `params_inf.m`, `baseline_*.m`, `run_powell.m`, `run_fmincon.m` | `draw_baseline_comparison.m` | Inputs and outputs under `baseline_comparison/results/`; expected outputs are loss curves, final-loss summaries, runtime/loss scatter, and trajectory/residual comparison panels |
| Normal-form examples | `Normal form/FHN/`, `Normal form/cancer/` | `run_fhn_fmam.m`, `run_cancer_fmam.m`, `build_*_normal_form_data.m` | `plot_*_normal_form.m`, `plot_*_fmam_amplitude_invariant.m`, `plot_*_fmam_period_invariant.m`, `reproduce_*_normal_form.m` | Inputs include checked-in `.mat` result files; expected outputs are normal-form and FMAM invariant panels |
| Repressilator supplementary panel | `RLT/` | `run_repressilator_fmam.m` | `draw_figS15.m` | Inputs include `figS15_repressilator_data.mat`; expected outputs are `Fig_S15.png` and `Fig_S15.pdf` |

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
