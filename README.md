
decorate: Differential Epigenetic Coregeulation Test

# Install
```r
library(devtools)

# first install sLED
install_github("lingxuez/sLED")

# then install decorate
install_github('https://github.com/GabrielHoffman/decorate.git')
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
sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute=20)

# get summary of results
df = summary( sledRes )
```
