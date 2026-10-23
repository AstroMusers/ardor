# _ardor_

Do stellar flares prefer a particular orbital phase, or are we seeing a pattern in a small sample? We use Ardor to detect and model flares in photometry, then test their timing for correlations. Our primary science case is magnetic star-planet interactions. We also look for links with stellar rotation, using simulations and statistical tests to examine how much the data can tell us.

## Repository paths

Set `ARDOR_PATH` to the checkout directory when using Ardor from another location. If it is unset, the package uses its checkout root. Local inputs and generated figures belong in the ignored `data/` and `visuals/` directories.

# Flare Detection Pipeline
`ardor` incorporates a flare detection pipeline to identify, characterize, and validate stellar flares in photometry. This pipeline uses four **tiers**:
- Tier 0: Light Curve Detrending
- Tier 1: Outlier identification
- Tier 2: Preliminary Model Parameterization
- Tier 3: Bayesian modeling using `allesfitter` or transdimensional inference with PCAT

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
renderer with `python -m pip install -e '.[bokeh]'`; the notebook then writes
hoverable HTML views of the flux posterior and flare-count posterior beside the
static figures under `examples/visuals/`. Its 1,200-sweep smoke chain is for
workflow illustration, not convergence or recovery-rate assessment.

## Testing flare recovery
We inject simulated flares into real light curves and run injection-recovery and precision-recall tests across flare parameter space. These tests help us see which flares the pipeline finds and which it misses.

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
Ardor supports the Kolmogorov-Smirnov (KS), Anderson-Darling (AD), and Kuiper (KU) tests, using `scipy`, `skgof`, and `astropy`, respectively. Each compares the empirical cumulative distribution function (eCDF) of the phase-folded flare epochs with the assumed uniform distribution. The resulting statistic has its own null distribution. Its p-value answers the question:

*If the flare phases were drawn from the assumed uniform distribution, how often would this test return a statistic at least as extreme as the one observed?*

That probability describes how unusual the statistic is under the uniform-phase assumption. It does not give the probability that the assumption is true. These tests assess departures from uniformity without fitting the shape of a phase-dependent flare distribution.

## Unbinned Likelihood Analysis
To examine the shape of the phase distribution, we fit the flare phases directly with a von Mises probability density function (PDF):

$$ f(x,\kappa) = \frac{\exp{\kappa \cos{x}}}{2\pi I_{0}(\kappa)}.$$
The von Mises PDF plays the role of a Gaussian on a circle. A cluster near the phase boundary stays one cluster instead of being split between the ends of a histogram. The location parameter shifts the peak around the circle, and the concentration parameter controls its width. Ardor fits these parameters by minimizing the negative log-likelihood with `scipy.optimize.minimize`.

![](../ardor/graphics/CDF+VM_Poster.png)
# Acknowledgements

The development of ardor has been supported by the McDonnell Center for Space Sciences at Washington University in St. Louis.
