

#' Extract sample-level correlation scores
#' 
#' Extract sample-level correlation scores for each cluster
#' 
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param gRanges GenomciRanges corresponding to the rows of epiSignal
#' @param clustList list of cluster assignments
#' @param method "deltaSLE", "Delaneau"
#' @param method.corr Specify type of correlation: "pearson", "kendall", "spearman"
#' @param BPPARAM parameters for parallel evaluation
#' 
#' @return matrix of scores of each sample for each cluster
#'
#' @examples
#' library(GenomicRanges)
#' 
#' # load data
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' # adjacentCount is the number of adjacent peaks considered in correlation
#' treeList = runOrderedClusteringGenome( simData, simLocation)
#' 
#' # Choose cutoffs and return cluster
#' treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c( 10, 20, 30, 40, 50) )
#' 
#' # Evaluate strength of correlation for each cluster
#' clstScore = scoreClusters(treeList, treeListClusters )
#' 
#' # Filter to retain only strong clusters
#' clustInclude = retainClusters( clstScore, "LEF", 0.30 )
#' 
#' # get retained clusters
#' treeListClusters_filter = filterClusters( treeListClusters, clustInclude)
#' 
#' # collapse similar clusters
#' treeListClusters_collapse = collapseClusters( treeListClusters_filter, simLocation)
#' 
#' # get correlation scores for each sample for each cluster
#' corScores = extractCorrelationScores( simData, simLocation, treeListClusters_collapse )
#' 
#' @seealso sle.score delaneau.score
#' @importFrom BiocParallel bpparam bplapply
#' @export
extractCorrelationScores = function(epiSignal, gRanges, clustList, method=c("deltaSLE", "Delaneau"), method.corr=c("pearson", "kendall", "spearman"), BPPARAM = bpparam()){

	method = match.arg( method )
	method.cor = match.arg( method.corr )

	corScores = lapply(names(clustList), function(id){
		res = lapply( names(clustList[[id]]), function(chrom){
			res = bplapply( unique(clustList[[id]][[chrom]]), function(idx){

				peakIDs = names(which(clustList[[id]][[chrom]] == idx))

				if( length(peakIDs) == 1){
					score = rep(NA, ncol(epiSignal))
					names(score) = colnames(epiSignal)
				}else if(  method == "deltaSLE"){
					score = sle.score( t(epiSignal[peakIDs,,drop=FALSE]), method.corr )
				}else if( method == "Delaneau" ){					
					score = delaneau.score( t(epiSignal[peakIDs,,drop=FALSE]), method.corr )
				}
				score
			}, BPPARAM=BPPARAM)
			res = do.call("rbind", res)

			# name rows by cluster id
			rownames(res) = paste(id, chrom, unique(clustList[[id]][[chrom]]), sep='_')

			# filter out rows with all NA's
			keep = apply(res, 1, function(x) any(!is.na(x)))
			res[keep,,drop=FALSE]
		})
		do.call("rbind", res)
	})
	
	# return matrix of clutser scores
	do.call("rbind", corScores)
}

#' Get genome coordinates for each cluster
#' 
#' Get genome coordinates for each cluster as a GRanges object
#' 
#' @param gRanges GenomciRanges corresponding to the rows of epiSignal
#' @param clustList list of cluster assignments
#' 
#' @return GRanges object
#' 
#' @importFrom GenomicRanges ranges
#'
#' @examples
#' library(GenomicRanges)
#' 
#' # load data
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' # adjacentCount is the number of adjacent peaks considered in correlation
#' treeList = runOrderedClusteringGenome( simData, simLocation)
#' 
#' # Choose cutoffs and return cluster
#' treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c( 10, 20, 30, 40, 50) )
#' 
#' # Get start and end coordinates for each cluster
#' # cluster name is parameter:chrom:cluster
#' getClusterRanges( simLocation, treeListClusters)
#' 
#' @export
getClusterRanges = function(gRanges, clustList){

	res = lapply(names(clustList), function(id){
		res = lapply( names(clustList[[id]]), function(chrom){
			res = lapply( unique(clustList[[id]][[chrom]]), function(idx){

				peakIDs = names(which(clustList[[id]][[chrom]] == idx))

				gr = ranges(gRanges[peakIDs])
				names(gr) = paste(id, chrom, idx, sep='_')
				gr
			})
			do.call('c', res)
		})
		do.call('c', res)
	})

	do.call('c', res)
}


#' Get name of each cluster
#' 
#' Get name of each cluster as parameter:chrom:cluster
#' 
#' @param clustList list of cluster assignments
#' 
#' @return array of cluster names
#'
#' @examples
#' library(GenomicRanges)
#' 
#' # load data
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' # adjacentCount is the number of adjacent peaks considered in correlation
#' treeList = runOrderedClusteringGenome( simData, simLocation)
#' 
#' # Choose cutoffs and return cluster
#' treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c( 10, 20, 30, 40, 50) )
#' 
#' # Name of each cluster is parameter:chrom:cluster
#' getClusterNames( treeListClusters)
#'
#' @export 
getClusterNames = function(clustList){

	res = lapply(names(clustList), function(id){
		res = lapply( names(clustList[[id]]), function(chrom){
			res = lapply( unique(clustList[[id]][[chrom]]), function(idx){

				peakIDs = names(which(clustList[[id]][[chrom]] == idx))

				c(parameter =id, chrom=chrom, cluster=idx, name = paste(id, chrom, idx, sep='_'))
			})
			do.call('rbind', res)
		})
		do.call('rbind', res)
	})
	res = data.frame(do.call('rbind', res), stringsAsFactors=FALSE)

	res$parameter = as.numeric(res$parameter)
	res$cluster = as.numeric(res$cluster)

	res
}

# Extract clutsters based on index
# Used for iterator
getClusterSubset = function(clustList, idx){

	# get identifies for each cluster
	df = getClusterNames( clustList )

	# get subset of cluster identifiers
	df_sub = df[idx,]

	# initialize count
	count <- 0L

	res = lapply(names(clustList), function(id){
		res = lapply( names(clustList[[id]]), function(chrom){
			res = lapply( unique(clustList[[id]][[chrom]]), function(idx){

				count <<- count + 1	

				clustID = paste(id, chrom, idx, sep='_')

				res = c()
				if( any(clustID == df_sub$name) ){
					peakInfo = clustList[[id]][[chrom]]
					res = peakInfo[peakInfo %in% idx]
					# res = list(res)					
					# names(res) = chrom
					# res = new('epiclustDiscreteList', res)
					res
				}
				res
			})	
		# drop empty entries
			res = Filter(length, res)	
			unlist(res)	
		})

		names(res) = names(clustList[[id]])
		new('epiclustDiscreteList', res)
	})

	names(res) = names(clustList)
	res = new('epiclustDiscreteListContain',res)

	# include result of getClusterNames to avoid re-running
	attr(res, 'getClusterNames') = df_sub

	res
}



















