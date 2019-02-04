
decorate: Differential Epigenetic Coregeulation Test

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
	install_github('https://github.com/GabrielHoffman/decorate.git', 
		build_vignettes=TRUE, dependencies=TRUE,
	  	repos=BiocInstaller::biocinstallRepos())
}
```

# Run example analysis
```r
library(GenomicRanges)
library(decorate)

data('decorateData')

# Evaluate hierarchical clsutering
treeList = runOrderedClusteringGenome( simData, simLocation )

# Choose cutoffs and return clutsers
treeListClusters = createClusters( treeList )

# Plot correlations and clusters in region defind by query
query = GRanges('chr1', IRanges(0, 1000))
plotDecorate( treeList, treeListClusters, query)

 # Simulate variable to split dataset by
set.seed(1)
metadata = data.frame( Disease = factor(sample(0:1, ncol(simData), replace=TRUE)))

# Evaluate Differential Correlation between two subsets of data
# use between 100 and 1000 permutations
sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute=c(100, 1000))

# get summary of results
df = summary( sledRes )
```
