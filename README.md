
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
	install_github('GabrielHoffman/decorate', 
		build_vignettes=TRUE, dependencies=TRUE,
		repos=BiocManager::repositories())
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
npermute = c(100, 1000)
sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute)

# get summary of results
df = summary( sledRes )
```

## examine top result
```r
# extract peak ID's from most significant cluster
peakIDs = getFeaturesInCluster( treeListClusters, df$chrom[1], df$cluster[1])

# plot comparison of correlation matrices for peaks in peakIDs
#  where data is subset by metadata$Disease
main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
plotCompareCorr( simData, peakIDs, metadata$Disease) + ggtitle(main)
```