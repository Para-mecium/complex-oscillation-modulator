# Longevity Model: SDE Noise Selection

## Scope

This note explains how to choose a reasonable noise level when converting the Sir2-HAP4 longevity oscillator to an SDE model. The target is not to claim a unique "true" sigma, but to choose a baseline that is consistent with the model class and with authoritative literature on synthetic gene oscillators and colored extrinsic noise.

## Recommended Baseline

- Noise form: multiplicative OU noise on kinetic rates, i.e. `rate -> rate * (1 + z_i)`.
- Correlation time: `tau_noise = 1 h` as a practical baseline for a slow extrinsic fluctuation.
- Baseline intensity: `sigma = 0.1`.
- Recommended sweep: `sigma = 0.05, 0.1, 0.15`.
- Stress test only: `sigma = 0.2`.

If the OU process is written as

```text
tau_noise dz = - z dt + sigma dW_t,
```

then its stationary variance is `sigma^2 / (2 * tau_noise)`. Therefore, when `tau_noise = 1 h`,

- `sigma = 0.1` gives `std(z) = 0.1 / sqrt(2) = 0.0707`, i.e. about a 7% relative fluctuation in the perturbed rate.
- `sigma = 0.15` gives about 10.6% relative fluctuation.
- `sigma = 0.2` gives about 14.1% relative fluctuation.

For a synthetic gene oscillator, this puts the baseline in the "moderate extrinsic noise" regime rather than in an aggressive stress-test regime.

## Why `sigma = 0.1` Is a Reasonable Default

### 1. The underlying deterministic oscillator is the synthetic Sir2-HAP4 longevity circuit

The deterministic model is not generic; it comes from the engineered Sir2-HAP4 negative-feedback oscillator introduced to slow cellular aging. That paper establishes the biological circuit and the kinetic scale that the SDE is perturbing.

Reference:

- Zhou Z, Liu Y, Feng Y, et al. Engineering longevity: design of a synthetic gene oscillator to slow cellular aging. Science. 2023;380:376-381. DOI: [10.1126/science.ade3403](https://doi.org/10.1126/science.ade3403)

## 2. For synthetic gene oscillators, fitted extrinsic noise is typically around 10%-15%, not 30%-40%

An experimentally grounded study of a synthetic gene oscillator separated intrinsic and extrinsic variability and fitted the extrinsic-noise parameter to data. Their fitted value was `Gamma = 0.12`, while the in silico validation example used `Gamma = 0.15`. Although `Gamma` in that paper is not the same symbol as the OU `sigma` here, both quantify the strength of slow cell-to-cell / cell-cycle-scale parameter variability in a synthetic oscillator. This puts a `sigma = 0.1` OU baseline in the right quantitative ballpark.

Reference:

- Cao AV, Clarke RW, Heltberg ML, et al. Sources of variability in a synthetic gene oscillator. PLoS Computational Biology. 2016;12:e1004674. DOI: [10.1371/journal.pcbi.1004674](https://doi.org/10.1371/journal.pcbi.1004674)

Useful quantitative statements from that paper:

- the fitted extrinsic-noise level was `Gamma = 0.12`;
- a toy-oscillator validation used `Gamma = 0.15`;
- the extrinsic variations evolved on a slow timescale of about 5 cell generations.

These numbers argue for treating about 10% noise as a realistic baseline and about 15% as the upper end of a typical range.

## 3. OU / colored extrinsic noise is the biologically appropriate modeling choice

For gene-expression systems, extrinsic fluctuations are usually not white and instantaneous. They are colored, nonspecific, and often persist on timescales comparable to the cell cycle. That is exactly the situation where OU noise is a standard reduced model.

Reference:

- Shahrezaei V, Ollivier JF, Swain PS. Colored extrinsic fluctuations and stochastic gene expression. Molecular Systems Biology. 2008;4:196. DOI: [10.1038/msb.2008.31](https://doi.org/10.1038/msb.2008.31)

So, for this longevity model, using OU noise on synthesis/degradation rates is not just numerically convenient; it matches the standard interpretation of extrinsic biological variability.

## 4. Robust-oscillation literature supports testing moderate parameter variability rather than only extreme noise

In the oscillator-robustness literature, extrinsic noise is commonly modeled as parameter variability, and robustness is assessed by how period and amplitude survive that variability. This supports evaluating the longevity oscillator under moderate parameter noise first, then under stronger stress-test conditions.

Reference:

- Qiao L, Zhang ZB, Zhao W, Wei P, Zhang L. Network design principle for robust oscillatory behaviors with respect to biological noise. eLife. 2022;11:e76188. DOI: [10.7554/eLife.76188](https://doi.org/10.7554/eLife.76188)

## Practical Recommendation For This Repo

For the current Sir2-HAP4 SDE implementation, the most defensible default is:

```text
tau_noise = 1 h
sigma = 0.1
```

and the most informative sensitivity sweep is:

```text
sigma in {0.05, 0.1, 0.15}
```

If a stronger-noise figure is needed, use

```text
sigma = 0.2
```

but label it explicitly as a robustness / stress-test case rather than as the baseline biological setting.

## Bottom Line

- `sigma = 0.1` is the best baseline choice.
- `sigma = 0.05-0.15` is the literature-consistent typical range.
- `sigma = 0.2` is acceptable for robustness checks.
- There is no strong literature basis for treating `sigma >= 0.3` as the default for this synthetic longevity oscillator.
