
# decorate: Differential Epigenetic Correlation Test
# decorate: _D_ifferential <u>E</u>pigenetic <u>Cor</u>rel<u>a</u>tion <u>Te</u>st
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

## [Vignette: run example analysis](http://deepfigv.mssm.edu/img/software/decorate/decorate_example.html)

## [Manual](http://deepfigv.mssm.edu/img/software/decorate/decorate-manual.pdf)

## [Simulations and properties of methods](http://deepfigv.mssm.edu/img/software/decorate/simulations.html)

## Data analysis from manuscript
 - [DNA methylation from kidney renal clear cell carcinoma](http://deepfigv.mssm.edu/img/software/decorate/KIRC.html)
 - [DNA methylation from schizophrenia brains](http://deepfigv.mssm.edu/img/software/decorate/methyl_scz_decorate.html)
 - [ATAC-seq from schizophrenia brains](http://deepfigv.mssm.edu/img/software/decorate/atac_local_corr.html)
 - [Histone modification ChIP-seq from human brain](http://deepfigv.mssm.edu/img/software/decorate/EpiMap.html)




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
