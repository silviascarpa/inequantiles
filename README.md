# inequantiles

<!-- badges: start -->
<!-- badges: end -->

**inequantiles** is an R package for estimating quantile-based inequality indicators from survey microdata, with full support for complex sampling designs and weighted quantile estimation.

📖 **Documentation and vignette**: [silviascarpa.github.io/inequantiles](https://silviascarpa.github.io/inequantiles/)

## Features

- **Quantile ratio index (QRI)**: estimation from survey data, parametric distributions, and grouped data
- **Traditional indicators**: quintile share ratio (QSR), Palma ratio, percentile ratios (e.g., P90/P10)
- **Weighted quantile estimation**: multiple interpolation rules (types 4–9 and Harrell-Davis) for complex survey data
- **Influence functions**: linearization-based techniques for QRI, QSR, Gini and quantiles
- **Grouped data**: quantiles, QRI and Gini from frequency tables (e.g., fiscal/administrative data)
- **Variance estimation**: rescaled bootstrap for complex sampling designs

## Installation

```r
# install.packages("devtools")
devtools::install_github("silviascarpa/inequantiles")
```

## Quick Start

```r
library(inequantiles)
data(synthouse)

# Weighted quantiles
csquantile(y       = synthouse$eq_income,
           weights = synthouse$weight,
           probs   = c(0.25, 0.5, 0.75),
           type    = 6)

# Quantile ratio index
qri(y       = synthouse$eq_income,
    weights = synthouse$weight)

# All indicators at once
inequantiles(y          = synthouse$eq_income,
             weights    = synthouse$weight,
             indicators = "all")

# Inequality curve
plot_inequality_curve(y       = synthouse$eq_income,
                      weights = synthouse$weight,
                      main    = "Inequality curve — synthouse")
```


## Variance Estimation

Standard errors for all indicators can be estimated simultaneously via the
rescaled bootstrap, using a single bootstrap loop for directly comparable results:

```r
inequantiles(
  y          = synthouse$eq_income,
  weights    = synthouse$weight,
  indicators = "all",
  se         = TRUE,
  data       = synthouse,
  strata     = "NUTS2",
  psu        = "municipality",
  B          = 200,
  seed       = 42
)
```

For custom estimators or more control over the bootstrap, use `rescaled_bootstrap()` directly.

## Influence Functions

Influence functions measure how much each observation affects an estimate —
useful for analytical variance estimation and diagnosing influential observations.

```r
# Influence function for the QRI
if_qri(y       = synthouse$eq_income,
       weights = synthouse$weight)

# Influence function for the QSR
if_qsr(y       = synthouse$eq_income,
       weights = synthouse$weight)

# Influence function for the Gini coefficient
if_gini(y       = synthouse$eq_income,
        weights = synthouse$weight)

# Influence function for the median
if_quantile(y       = synthouse$eq_income,
            weights = synthouse$weight,
            probs   = 0.5)
```


## Grouped Data

When only frequency tables are available (e.g., tax records):

```r
income_freq  <- c(120, 180, 150, 80, 40, 20, 10)
income_tot   <- c(18800, 16300, 44700, 33900, 21500, 22100, 98300)
income_lower <- c(0, 15000, 30000, 45000, 60000, 80000, 100000)
income_upper <- c(15000, 30000, 45000, 60000, 80000, 100000, 150000)

quantile_grouped(freq = income_freq,
                 lower_bounds = income_lower,
                 upper_bounds = income_upper,
                 probs = c(0.25, 0.5, 0.75))

qri_grouped(freq = income_freq,
            lower_bounds = income_lower,
            upper_bounds = income_upper)

gini_grouped(Y = income_tot, freq = income_freq)
```


## Citation

If you use **inequantiles** in your research, please cite:

> Scarpa, S. (2025). *inequantiles: Quantile-Based Inequality Measures for Survey Data*. R package. https://github.com/silviascarpa/inequantiles

## Getting Help

- 📋 [Open an issue](https://github.com/silviascarpa/inequantiles/issues)
- 📧 silvia.scarpa@unimore.it
