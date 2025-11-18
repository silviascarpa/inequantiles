# inequantiles

<!-- badges: start -->
<!-- badges: end -->

**inequantiles** provides R functions for the estimation of quantile-based inequality indicators, 
such as the *quantile ratio index (QRI)*, and quantiles, also allowing for complex sampling data.

## Installation

You can install the development version of **inequantiles** from GitHub with:

```r
# install.packages("devtools")   # run this line if you don't have it yet
devtools::install_github("silviascarpa/inequantiles")
```

## Example

```r
library(inequantiles)

# Example: compute a quantile ratio index (QRI)
set.seed(123)
N <- 1000
y <- rlnorm(N)  # simulated income data

qri(y) ## Estimate the QRI on a srs

# On (synthetic) survey data

data(synthouse) # Import household survey dataset

csquantile(synthouse$eq_income, weight = synthouse$weight, probs = c(0.25, 0.50, 0.75), type = 6) ### Estimate quantiles on survey data

## QRI estimation on survey data with type 6
qri(y = synthouse$eq_income, weight = synthouse$weight, 
type = 4) 

## Palma index estimation on survey data
palma_ratio(y = synthouse$eq_income, weight = synthouse$weight, 
type = 6) 

## P90/P10 ratio estimation on survey data
ratio_quantiles(y = synthouse$eq_income, weight = synthouse$weight, 
prob_numerator = 0.90, prob_denominator = 0.10) 

## Influence function estimation of the QRI (on some observations, as example)
if_qri(y = synthouse$eq_income[1:100], weight = synthouse$weight[1:100])

```

