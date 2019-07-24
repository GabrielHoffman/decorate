
# decorate: Differential Epigenetic Coregeulation Test

![alt text](https://hoffmg01.u.hpc.mssm.edu/software/decorate/workflow.png)

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

# Dependencies
Depending on your system, you may need to install these dependencies first: 
- udunits 
- proj
- gdal 
- geos

###### Mac OS X
```brew install udunits proj gdal geos```

###### Ubuntu
```sudo apt-get install libudunits2-dev libproj-dev proj-data proj-bin libgeos-dev libgeos-c1v5 libgdal-dev```

###### Windows
I have not tried to install these dependencies on Windows, but decorate should install on Windows with no issues.


## Vignette: [Run example analysis](https://hoffmg01.u.hpc.mssm.edu/software/decorate/decorate_example.html)
