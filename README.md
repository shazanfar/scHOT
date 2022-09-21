# scHOT: single-cell higher order testing

<img src="man/figures/scHOT_hex.png" align="right" width="250"/>

## Installation

`scHOT` is available as a [Bioconductor package](https://bioconductor.org/packages/release/bioc/html/scHOT.html)

For the latest version, install the following packages using `BiocManager`:

```r
# install.packages("BiocManager")
BiocManager::install(c("S4Vectors", "SummarizedExperiment", "SingleCellExperiment", 
"Matrix", "IRanges", "BiocParallel", "reshape", "ggplot2", "igraph", "grDevices", "ggforce"))
```

Then install the latest version of `scHOT` using `devtools`:

```r
library(devtools)
devtools::install_github("shazanfar/scHOT")
```

## Vignette

You can find the latest vignette at this website: https://shazanfar.github.io/scHOT/.


## Contact

shila.ghazanfar \<at\> sydney.edu.au


