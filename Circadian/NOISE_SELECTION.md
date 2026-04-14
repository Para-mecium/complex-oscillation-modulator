# Circadian Model: SDE Noise Selection

## Scope

This note explains how to choose a reasonable noise level when converting the current circadian model to an SDE. The model in this repository is a minimal Kim-Forger-type clock, so the most relevant quantitative evidence comes from stochastic Kim-Forger-based circadian studies rather than from unrelated oscillator classes.

## Recommended Baseline

- Noise form in this repo: multiplicative OU noise on the three dynamic rates.
- Correlation time: `tau_noise = 1 h`.
- Baseline intensity: `sigma = 0.1`.
- Recommended sweep: `sigma = 0.05, 0.1, 0.15`.
- Strong-noise regime: `sigma = 0.2-0.3`.
- Stress test only: `sigma = 0.4`.

If the OU process is written as

```text
tau_noise dxi = - xi dt + sigma dW_t,
```

then its stationary variance is `sigma^2 / (2 * tau_noise)`. Therefore, when `tau_noise = 1 h`,

- `sigma = 0.1` gives `std(xi) = 0.0707`, i.e. about a 7% relative fluctuation in each perturbed rate.
- `sigma = 0.2` gives about a 14% relative fluctuation.
- `sigma = 0.4` gives about a 28% relative fluctuation.

For a minimal circadian clock, `sigma = 0.4` is therefore not a mild baseline. It is a deliberately strong perturbation.

## Why `sigma = 0.1` Is the Best Baseline

### 1. The deterministic core is a Kim-Forger circadian oscillator

The current circadian model belongs to the Kim-Forger family based on protein sequestration, which is the right reference class for quantitative noise selection.

Reference:

- Kim JK, Forger DB. A mechanism for robust circadian timekeeping via stoichiometric balance. Molecular Systems Biology. 2012;8:630. DOI: [10.1038/msb.2012.62](https://doi.org/10.1038/msb.2012.62)

## 2. A stochastic Kim-Forger-based circadian model fitted to data gives a control noise intensity of about 0.1

A later stochastic Kim-Forger-based study fitted the noise strength directly to zebrafish circadian reporter data. In the control condition, the fitted value was

```text
sigma = 0.098 +- 0.006
```

and representative treatment conditions gave values such as

- `sigma = 0.124 +- 0.016` for a noisier condition;
- `sigma = 0.023 +- 0.001` for a quieter condition.

This is the strongest quantitative evidence available for a minimal Kim-Forger-type stochastic clock. It says that a baseline around `0.1` is typical, while `0.4` is much larger than the fitted control level.

Reference:

- Kumpošt V, Bera K, Manjunath H, et al. A stochastic oscillator model simulates the entrainment of vertebrate cellular clocks by light. Scientific Reports. 2021;11:13297. DOI: [10.1038/s41598-021-93913-2](https://doi.org/10.1038/s41598-021-93913-2)

Important interpretation point:

- that paper uses additive Wiener noise rather than multiplicative OU noise;
- however, it is still the closest quantitative prior because it is built from the same minimal Kim-Forger clock family and is fitted to data.

Therefore, for this repository's multiplicative OU formulation, `sigma = 0.1` is the most defensible baseline translation.

## 3. In stochastic Kim-Forger entrainment studies, too much noise is not uniformly beneficial

Another Kim-Forger-based stochastic study showed that there is an optimal noise intensity for entrainment: increasing noise can help up to a point, but excessive noise degrades entrainment again. This supports a workflow where one chooses a moderate baseline and then explores sensitivity, instead of setting an overly large default.

Reference:

- Kumpošt V, Hilbert L, Mikut R. Noise facilitates entrainment of a population of uncoupled limit cycle oscillators. Journal of the Royal Society Interface. 2023;20:20220781. DOI: [10.1098/rsif.2022.0781](https://doi.org/10.1098/rsif.2022.0781)

## 4. Colored / correlated extrinsic noise remains a biologically sensible choice

Even though the best-fitted Kim-Forger circadian paper above uses additive Wiener noise, correlated extrinsic variability is still a standard and biologically meaningful reduced model for gene-regulatory systems. If the goal here is robustness analysis under parameter fluctuations, OU noise is appropriate.

Reference:

- Shahrezaei V, Ollivier JF, Swain PS. Colored extrinsic fluctuations and stochastic gene expression. Molecular Systems Biology. 2008;4:196. DOI: [10.1038/msb.2008.31](https://doi.org/10.1038/msb.2008.31)

## Practical Recommendation For This Repo

For the current circadian SDE implementation, the most defensible default is:

```text
tau_noise = 1 h
sigma = 0.1
```

and the most informative sensitivity sweep is:

```text
sigma in {0.05, 0.1, 0.15}
```

If a stronger robustness figure is desired, use

```text
sigma = 0.2 or 0.3
```

If `sigma = 0.4` is kept, it should be described as

```text
a strong-noise robustness test
```

not as the default "typical" biological setting.

## Relation To The Current Repository

Because the fitted control value in the closest Kim-Forger stochastic literature is about `0.1`, the repository should be internally interpreted as follows:

- a default of `sigma = 0.1` is literature-aligned;
- a choice of `sigma = 0.4` is substantially stronger and should be labeled accordingly.

## Bottom Line

- `sigma = 0.1` is the best baseline choice.
- `sigma = 0.05-0.15` is the most defensible typical range.
- `sigma = 0.2-0.3` is a strong-noise regime.
- `sigma = 0.4` should be treated as a stress test, not as a typical baseline.
