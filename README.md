# methyLImp2
**methyLImp2** package is the implementation of the updated version of the _methyLImp_ method for missing value imputation in methylation data. Compared to the previous version, implementation is parallelised over chromosomes since probes on different chromosomes are usually independent. Moreover, mini-batch approach to reduce the runtime in case of large number of samples is available.

## Installation

Package is available from [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/methyLImp2.html)! To install, run
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methyLImp2")
```

To install from Github, run
```
#install.packages("devtools")

library("devtools")
install_github("annaplaksienko/methyLImp2", build_vignettes = TRUE)
```

## Vignette
To open the vignette, run

```{r run methyLImp}
library(methyLImp2)
browseVignettes("methyLImp2")
```

