
# decorate: Differential Epigenetic Coregeulation Test

![alt text](https://hoffmg01.u.hpc.mssm.edu/software/decorate/workflow.png)

# Dependencies
Depending on your system, you may need to install dependencies first. Specifically, install:
- udunits 
- proj
- gdal 
- geos

###### Mac OS X
```brew install udunits proj gdal geos```

###### Ubuntu
```sudo apt-get install libproj-dev proj-data proj-bin libgeos-dev libgeos-c1v5 libgdal-dev```

I have not tried to install these dependencies on Windows

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