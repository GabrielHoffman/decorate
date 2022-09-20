
# *D*ifferential *E*pigenetic *Cor*rel*a*tion *Te*st
![](https://hoffmg01.u.hpc.mssm.edu/software/decorate/workflow.png)

# Install
```r
library(devtools)

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
Depending on your system, you may need to install these system dependencies first: 

- udunits 
- proj
- gdal 
- geos

###### Mac OS X
```
brew install udunits proj gdal geos
```

###### Ubuntu
```
apt-get install libudunits2-dev libproj-dev proj-data proj-bin libgeos-dev libgeos-c1v5 libgdal-dev
```

###### CentOS
```
yum install udunits2-devel
yum install proj-devel
yum install gdal-devel
yum install geos-devel
```

###### Windows
decorate should install on Windows without needing to install these dependencies.


### Data analysis from manuscript
 - [DNA methylation from kidney renal clear cell carcinoma](https://ghoffman-cdn.s3.us-east-2.amazonaws.com/decorate/KIRC.html)
 - [DNA methylation from schizophrenia brains](https://ghoffman-cdn.s3.us-east-2.amazonaws.com/decorate/methyl_scz_decorate.html)
 - [ATAC-seq from schizophrenia brains](https://ghoffman-cdn.s3.us-east-2.amazonaws.com/decorate/atac_local_corr.html)
 - [Histone modification ChIP-seq from human brain](https://ghoffman-cdn.s3.us-east-2.amazonaws.com/decorate/EpiMap.html)

### Simulations and properties of methods
  - [Simulations](https://ghoffman-cdn.s3.us-east-2.amazonaws.com/decorate/simulations.html)



<!--
## [Vignette: run example analysis](https://hoffmg01.u.hpc.mssm.edu/software/decorate/decorate_example.html)

## [Manual](https://hoffmg01.u.hpc.mssm.edu/software/decorate/decorate-manual.pdf)

## [Simulations and properties of methods](https://hoffmg01.u.hpc.mssm.edu/software/decorate/decorate_example.html)

## Data analysis from manuscript
 - [DNA methylation from kidney renal clear cell carcinoma](https://hoffmg01.u.hpc.mssm.edu/software/decorate/KIRC.html)
 - [DNA methylation from schizophrenia brains](https://hoffmg01.u.hpc.mssm.edu/software/decorate/methyl_scz_decorate.html)
 - [ATAC-seq from schizophrenia brains](https://hoffmg01.u.hpc.mssm.edu/software/decorate/atac_local_corr.html)
 - [Histone modification ChIP-seq from human brain](https://hoffmg01.u.hpc.mssm.edu/software/decorate/EpiMap.html)
-->

  [Results files](https://www.synapse.org/#!Synapse:syn20742092)

  [Analysis code](https://github.com/GabrielHoffman/decorate_analysis)
