# changepointsR

Changepoint algorithms and inference for change in mean model. Paper: https://arxiv.org/pdf/2301.05636.pdf

********************************************************************************************************************************************
## Installation

To install this package:
```
devtools::install_github("rachelcarrington/changepointsR")
```

If you want to use L0 inference, you will also need to download the package ChangepointInference (see https://github.com/jewellsean/ChangepointInference):
```
devtools::install_github("jewellsean/ChangepointInference")
```

To use seeded binary segmentation, you will need to download the function `seedBS` from https://github.com/kovacssolt/SeedBinSeg.

Other dependencies (CRAN packages):
* ggplot2 (for plots)
* changepoint (for data)

********************************************************************************************************************************************

### Changepoint algorithms
The following changepoint algorithms are included:
* binary segmentation: `binary_segmentation`
* wild binary segmentation: `wild_binary_segmentation`
* seeded binary segmentation: `wild_binary_segmentation` with `seeded = TRUE`
* narrowest over threshold: `narrowest_over_threshold`

********************************************************************************************************************************************

### Post-selection inference
To do post-selection inference:
* Use `calculate_pvals` for binary segmentation, wild/seeded binary segmentation, and narrowest over threshold. </br>
* Use `l0_segmentation_psi` for L0 segmentation.

********************************************************************************************************************************************

## Examples

Binary segmentation:

```
### Generate some data
set.seed(100)
y <- rnorm(200) + c(rep(1,40), rep(-1,40), rep(1,40), rep(-1,40), rep(1,40))

### Implement binary segmentation
results <- binary_segmentation(y, threshold=4)

### Calculate p-values
pvals <- calculate_pvals(y, method="bs", results=results, N=10, h=10, return_probs=TRUE)
print(pvals)
```

L0 segmentation:

```
### Generate some data
set.seed(100)
y <- rnorm(200) + c(rep(1,40), rep(-1,40), rep(1,40), rep(-1,40), rep(1,40))

### Implement L0 segmentation
results <- changepoint_estimates(y, "L0", 10)

### Calculate p-values
pvals <- l0_segmentation_psi(y, lambda=10, N=10, h=10, sigma2=1, include_original=TRUE, num_pvals=1)
print(pvals)
```
