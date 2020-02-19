# scHOT: single-cell higher order testing



## Installation


Install the following packages using `BiocManager`:

```r
# install.packages("BiocManager")
BiocManager::install(c("S4Vectors", "SummarizedExperiment", "SingleCellExperiment", 
"Matrix", "IRanges", "BiocParallel", "reshape", "ggplot2", "igraph", "grDevices", "ggforce"))
```

    
    
Then install `scHOT` using `devtools`:

```r
library(devtools)
devtools::install_github("shazanfar/scHOT")
```


## Vignette

You can find the vignette at this website: https://shazanfar.github.io/scHOT/.





