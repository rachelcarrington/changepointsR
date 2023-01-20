# changepointsR

Changepoint algorithms and inference for change in mean model. Paper: https://arxiv.org/pdf/2301.05636.pdf

The following changepoint algorithms are included:
* binary segmentation (binary_segmentation)
* wild binary segmentation

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

To use seeded binary segmentation, you will need to download the function \code{seedBS}


********************************************************************************************************************************************



Changepoint algorithms: </br>
binary_segmentation </br>
wild_binary_segmentation </br>
narrowest_over_threshold

Post-selection inference: </br>
Use calculate_pvals for binary segmentation, wild binary segmentation, and narrowest over threshold. </br>
Use l0_segmentation_psi for L0 segmentation.

Dependencies: </br>
ChangepointInference (https://github.com/jewellsean/ChangepointInference) for L0 segmentation

For plots: </br>
ggplot2 </br>
Changepoint (for data)
