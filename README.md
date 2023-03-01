# methyLImp2
**methyLImp2** package is the implementation of the updated version of the _methyLImp_ method for missing value imputation in methylation data. Compared to the previous version, implementation is parallelised over chromosomes since probes on different chromosomes are usually independent. Moreover, mini-batch approach to reduce the runtime in case of large number of samples is available.

## Installation
To install from github, run
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

