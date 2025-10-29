# _ardor_

**Welcome to ardor**

`ardor` is a package dedicated to detecting, modeling, and analyzing correlated flares in photometry. It featues a multi-tiered flare detection pipeline, simulation tools, statistical testing, and forward modeling functions, as well as a visualization suite. While its primary goal is to detect correlated flare arising from magnetic star-planet interactions, `ardor` also searches for correlation with stellar rotation.

# Flare Detection Pipeline
`ardor` incorporates a flare detection pipeline to identify, characterize, and validate stellar flares in photometry. This pipeline uses four **tiers**:
- Tier 0: Light Curve Detrending
- Tier 1: Outlier identification
- Tier 2: Preliminary Model Parameterization
- Tier 3: Bayesian Modeling using `allesfitter`

The end product of the `ardor` pipeline is a list of identified flares and flare parameters, including epoch, flare ampltidue, duration, and bolometric energy. This pipeline can be found in the `Flare.py` module.

![](../ardor/graphics/Ardor_T1_T2_Rescale.jpg)

![](../ardor/graphics/Allesfitter.png)

## Pipeline Diagnostics
`ardor` supports injection-recovery and precision-recall tests using simulated flares in real light curves to check performance over flare paramter space.

![](../ardor/graphics/Injection_Recovery_Parameters.png)

# Statistical Tests
`ardor` can anaylze flare samples for a particular planetary system and determine if the flare samples are consistent with significant flare clustering. This is performed using two types of tests:
- Goodness-of-fit (GoF) tests
- Unbinned Likelihood Analysis

Both take as input phase folded flare epochs with the planetary orbital period and output their relevant test statistic and significance metric.

## GoF Tests
Three different GoF tests are supported: the Kolmolgorov-Smirnov (KS) test, the Anderson-Darling (AD) test, and the Kuiper (KU) test. These statistical tests are wrapped from the `scipy`, `skgof`, and `astropy` packages, respectively. `ardor` constructs the empirical cumulative density function (eCDF) from the given phase-folded flare epochs. The GoF tests compare the eCDF to the expected uniform distribution and determines a distance statistic between the eCDF and the uniform CDF. Each distance statistic is drawn from its own distribution which returns a p-value, which can be used to determine the significance of the flare sample. The p-value answers the question:

*What is the probability that the empirical sample is drawn from a uniform distribution?*

Thus, the GoF tests can answer if a sample is non-uniform to some signifiance level, but does not inform the shape of the sample.

## Unbinned Likelihood Analysis
To compliment the GoF tests, `ardor` uses a method called Unbinned Likelihood Analysis which aims to find the best-fit parameters that minimizes the likelihood function. Here, we directly fit a probability density function (PDF), the von Mises PDF:

$$ f(x,\kappa) = \frac{\exp{\kappa \cos{x}}}{2\pi I_{0}(\kappa)}.$$
The von Mises PDF can be thought of as a periodic Gaussian; it is unimodal and can cross periodic boundaries, ideal characterstics to find flare clustering a particular orbital phase. The likelihood function is minimized using `scipy.optimize.minimize`.

![](../ardor/graphics/CDF+VM_Poster.png)
# Acknowledgements

The development of ardor has been supported by the McDonnell Center for Space Sciences at Washington University in St. Louis.
