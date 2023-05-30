# scHOT: single-cell higher order testing

<img src="man/figures/scHOT_hex.png" align="right" width="250"/>

## Installation

`scHOT` is available as a [Bioconductor package](https://bioconductor.org/packages/release/bioc/html/scHOT.html)

For the latest release, install `scHOT` using `BiocManager`:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scHOT")
```

Alternatively, install the development version of `scHOT` using `devtools`:

```r
library(devtools)
devtools::install_github("shazanfar/scHOT")
```

## Vignette

You can find the latest vignette at this website: https://shazanfar.github.io/scHOT/.


## Contact

shila.ghazanfar \<at\> sydney.edu.au

## Citation

Ghazanfar, S., Lin, Y., Su, X. *et al.* Investigating higher-order interactions in single-cell data with scHOT. *Nat Methods* **17**, 799--806 (2020). https://doi.org/10.1038/s41592-020-0885-x
