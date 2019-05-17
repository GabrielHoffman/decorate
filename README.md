
# decorate: Differential Epigenetic Coregeulation Test

![alt text](https://hoffmg01.u.hpc.mssm.edu/software/decorate/workflow.png.png)

# Install
```r
library(devtools)

# first install sLED
install_github("lingxuez/sLED")

# Install decorate
# 	first, check for Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)){
	cat("Please install Bioconductor before continuing:\n")
	cat("see http://bioconductor.org/install/\n\n")
}else{
	install_github('GabrielHoffman/decorate', repos=BiocManager::repositories())
}
```

## Vignette: [Run example analysis](https://hoffmg01.u.hpc.mssm.edu/software/decorate/decorate_example.html)