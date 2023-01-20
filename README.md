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

## Example

```{r}
y = rnorm(100)
plot(y)
```
