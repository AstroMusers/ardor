# _ardor_

Do stellar flares favor a particular orbital phase, or can a small sample make coincidence look like a pattern? Ardor detects and models flares in photometry so we can test their timing against planetary orbits and stellar rotation. Our main science case is magnetic star-planet interactions. Simulations show how readily small samples can mimic a preferred phase.

## Repository paths

Set `ARDOR_PATH` to the checkout directory when using Ardor from another location. If it is unset, the package uses its checkout root. Local inputs and generated figures belong in the ignored `data/` and `visuals/` directories.

# Flare Detection Pipeline
`ardor` incorporates a flare detection pipeline to identify, characterize, and validate stellar flares in photometry. The pipeline has four tiers.
- **Tier 0.** Light-curve detrending
- **Tier 1.** Outlier identification
- **Tier 2.** Preliminary model parameterization
- **Tier 3.** Bayesian modeling with `allesfitter` or transdimensional inference with PCAT

The pipeline returns identified flares and their parameters, including epoch, amplitude, duration, and bolometric energy. The implementation is in the `Flare.py` module.

![](../ardor/graphics/Ardor_T1_T2_Rescale.jpg)

![](../ardor/graphics/Allesfitter.png)

## Variable-count PCAT fits

The existing `tier3()` pipeline uses Allesfitter to compare fixed models with up to three flare peaks. For a variable number of flares, install the optional PCAT backend and run the Tier 2 cutouts through `tier3_pcat()`:

```bash
pip install -e '.[pcat]'
```

```python
from ardor.Flares.pcat_flares import tier3_pcat

fits = tier3_pcat(
	tier_2_output_dir="path/to/tier2/flare_csvs",
	tier_3_output_dir="path/to/pcat/results",
	host_name="TIC_123456",
	numbswep=20000,
	max_flares=10,
)
```

Each flare cutout gets its own PCAT run directory. The PCAT posterior retains the sampled flare count and catalog samples instead of forcing them into the fixed-count Allesfitter CSV schema. A simulated end-to-end example is available with `python examples/pcat_variable_flares.py --smoke`.

The [variable-flare notebook](examples/pcat_variable_flares.ipynb) runs the
same synthetic Tier 2 input through `tier3_pcat()` and displays the full-baseline
fit and sampled flare-count distribution. Install the optional interactive
renderer with `python -m pip install -e '.[bokeh]'`. The notebook then writes
hoverable HTML views of the flux posterior and flare-count posterior beside the
static figures under `examples/visuals/`. The notebook runs 1,200 sweeps.
Convergence and recovery rates require separate analyses.

## Testing flare recovery
We inject simulated flares into real light curves and measure recovery and precision-recall across flare parameter space. These results show which signals the pipeline detects and which it misses.

![](../ardor/graphics/Injection_Recovery_Parameters.png)

## Roman flare fits

Install the Roman sampling dependencies with `pip install -e '.[roman]'`.
The MCMC and nested samplers share the flare model and CSV loader in
`ardor.Roman.flare_model`. Both build custom KDE priors with
`tdpy.numerics.build_kde_prior_from_samples`, retaining their existing grid
normalization and prior metadata.

# Statistical Tests
For a particular planetary system, we test whether the observed flare phases are consistent with a uniform distribution. Ardor provides two approaches:
- Goodness-of-fit (GoF) tests
- Unbinned Likelihood Analysis

Both use flare epochs folded on the planetary orbital period and return test statistics. Interpreting them also requires considering observing coverage and the assumptions of the chosen test.

## GoF Tests
Ardor supports the Kolmogorov-Smirnov (KS), Anderson-Darling (AD), and Kuiper (KU) tests, using `scipy`, `skgof`, and `astropy`, respectively. Each compares the empirical cumulative distribution function (eCDF) of the phase-folded flare epochs with the assumed uniform distribution. The statistic has a corresponding null distribution. Under the uniform-phase assumption, the p-value is the probability of obtaining a statistic at least as extreme as the observed one. A small p-value is evidence against uniform timing. To describe a phase-dependent pattern, fit the flare phases directly.

## Unbinned Likelihood Analysis
We fit the phase distribution with a von Mises probability density function (PDF):

$$ f(x \mid \mu, \kappa) = \frac{\exp\left[\kappa \cos(x-\mu)\right]}{2\pi I_{0}(\kappa)}. $$
The von Mises distribution is the circular counterpart of a Gaussian. Phases just below 1 and just above 0 are neighbors on the orbit, so flares clustered across that boundary remain together. The location parameter sets the peak's phase, and the concentration parameter sets how tightly the flares cluster around it. Ardor finds these parameters by minimizing the negative log-likelihood with `scipy.optimize.minimize`.

![](../ardor/graphics/CDF+VM_Poster.png)
# Acknowledgements

The development of ardor has been supported by the McDonnell Center for Space Sciences at Washington University in St. Louis.
