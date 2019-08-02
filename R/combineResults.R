


#' Combine results into a single data.frame
#'
#' Combine results into a single data.frame for easy post processing
#'
#' @param sledRes sLEDresults from evalDiffCorr()
#' @param clstScore cluster summary statistics from from scoreClusters()
#' @param treeListClusters epiclustDiscreteListContain from createClusters()
#' @param peakLocations GenomeRanges object
#' @param verbose show messages
#'
#' @examples
#' library(GenomicRanges)
#' library(EnsDb.Hsapiens.v86)
#'
#' # load data
#' data('decorateData')
#'
#' # load gene locations
#' ensdb = EnsDb.Hsapiens.v86
#' 
#' # Evaluate hierarchical clsutering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c( 10, 20) )
#' 
#' # Evaluate strength of correlation for each cluster
#' clstScore = scoreClusters(treeList, treeListClusters )
#' 
#' # Filter to retain only strong clusters
#' # If lead eigen value fraction (LEF) > 30% then keep clusters
#' # LEF is the fraction of variance explained by the first eigen-value
#' clustInclude = retainClusters( clstScore, "LEF", 0.30 )
#' 
#' # get retained clusters
#' treeListClusters_filter = filterClusters( treeListClusters, clustInclude)
#' 
#' # collapse redundant clusters
#' treeListClusters_collapse = collapseClusters( treeListClusters_filter, simLocation, jaccardCutoff=0.9)
#'
#' # Evaluate Differential Correlation between two subsets of data
#' sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters_collapse, npermute=c(20, 200, 2000))
#'
#' # Combine results for each cluster
#' df_results = combineResults( sledRes, clstScore, treeListClusters, simLocation)
#'
#' @export
combineResults = function( sledRes, clstScore, treeListClusters, peakLocations, verbose=TRUE){
  
  # get test results
  if( verbose ){
  	cat("Summarizing analysis...\n")
  }
  df = summary(sledRes)

  # combine with cluster summary statistics
  if( verbose ){
  	cat("Summarizing cluster properties...\n")
  }
  df_combine = merge( df, do.call('rbind', clstScore ), by=c("id", 'chrom', 'cluster') )

  # extract cluster locations
  if( verbose ){
	 cat("Collecting cluster locations...\n")
	}
  # res = lapply( seq_len(nrow(df_combine)), function(i){
  #   peakIDs = getFeaturesInCluster( treeListClusters, df_combine$chrom[i], df_combine$cluster[i], df_combine$id[i])
  #   data.frame( chrom = df_combine$chrom[i], 
  #               cluster = df_combine$cluster[i], 
  #               id = df_combine$id[i], 
  #               range(peakLocations[peakIDs]))
  # })
  # res = do.call("rbind", res)

  clustInfo = getClusterRanges( peakLocations, treeListClusters, verbose=verbose)
  clustDf = data.frame(clustInfo, stringsAsFactors=FALSE)
  colnames(clustDf)[colnames(clustDf)=="seqnames"] = "chrom"
  clustDf$chrom = as.character(clustDf$chrom ) 
  
  # merge
  if( verbose ){
	cat("Merging results...\n")
	}
  df_combine = merge( df_combine, clustDf[,c('chrom', 'cluster', 'id', 'start', 'end', 'width')], by=c("id", 'chrom', 'cluster') )

  # sort
  df_combine = df_combine[order(df_combine$pValue, -abs(df_combine$stat)),]
 
  rownames(df_combine) = c()

  df_combine
}




