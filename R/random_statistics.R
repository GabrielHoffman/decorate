


#' Run hierarchical clustering permuting features
#' 
#' Run hierarchical clustering permuting features to get statistics under the null
#' 
#' @param X data matrix were *rows* are features in sequential order
#' @param gr GenomicRanges object with entries corresponding to the *rows* of X
#' @param method 'adjclust': adjacency constrained clustering.  'hclustgeo': incorporate data correlation and distance in bp
#' @param quiet suppress messages
#' @param alpha use by 'hclustgeo': mixture parameter weighing correlations (alpha=0) versus chromosome distances (alpha=1)
#' @param adjacentCount used by 'adjclust': number of adjacent entries to compute correlation against
#' @param setNANtoZero replace NAN correlation values with a zero
#' @param method.corr Specify type of correlation: "pearson", "kendall", "spearman"
#' @param meanClusterSize select target mean cluster size.  Can be an array of values
#' 
#' @return list of clusterScores and cutoff values at 5% false positive rate
#'
#' @examples
#' library(GenomicRanges)
#'
#' # load data
#' data('decorateData')
#' 
#' # First, analysis of original data
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList, method='meanClusterSize', meanClusterSize=c(5, 10) )
#'
#' # Evaluate score for each cluster
#' clstScore = scoreClusters(treeList, treeListClusters )
#'
#' # Then, analysis of permuted data
#' # Evaluate hierarchical clustering
#' res = runPermutedData( simData, simLocation, meanClusterSize=c(5, 10)  ) 
#'
#' # LEF values for permuted data at 5% false positive rate
#' res$cutoffs$LEF
#'
#' # Retain clusters that pass this criteria
#' clustInclude = retainClusters( clstScore, "LEF", res$cutoffs$LEF )
#'
#' @export
runPermutedData = function( X, gr, method = c("adjclust", 'hclustgeo'), quiet=FALSE, alpha=0.5, adjacentCount=500, setNANtoZero=FALSE, method.corr = c("pearson", "spearman"), meanClusterSize=c(5, 10) ){

	# permute features (i.e. rows of dataset)
	idx = sample.int(nrow(X), nrow(X))

	X_perm = X[idx,]
	rownames(X_perm) = rownames(X)
	rm(X)

	# perform clustering and score each cluster
	treeList = runOrderedClusteringGenome( X_perm, gr, method, quiet, alpha, adjacentCount, setNANtoZero, method.corr)

	treeListClusters = createClusters( treeList, method='meanClusterSize', meanClusterSize )

	# get total number of clusters
	n_clusters = countClusters( treeListClusters )

	# score each cluster
	clstScore = scoreClusters(treeList, treeListClusters )

	# get cutoffs at 0.05 False Positive Rate
	cutoffs = lapply(clstScore, function(df){
		apply(df[df$N>1,c('mean_abs_corr', 'quantile75', 'quantile90', 'quantile95', 'LEF')], 2, function(x){
			x - x
			sort(x)[nrow(df)*0.95]
		})
		})
	cutoffs = do.call('rbind', cutoffs)

	# return score of each clutser
	list(clusterScores = clstScore, cutoffs = as.data.frame(cutoffs))
}



 

