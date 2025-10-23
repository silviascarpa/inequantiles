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
y <- rlnorm(100)  # simulated income data
w <- rlnorm(6, 0.4, length(y)) # simulated sampling weights
qri(y, w)

if_qri(y, w)
```
