#' Compute distance between peaks
#'
#' Given a query set of genome intervals, report distance to all others within window
#'
#' @param query GRanges object of intervals to query
#' @param windowSize with of window around each interval in bp
#' 
#' @details A smaller window size will create a covariance matrix that is faster to evaluate, but larger windows give a better approximation and are less likely to have negative eigen value
#'
#' @return for all pairs of peaks within windowSize, report distance 
#' 
#' @examples
#' library(GenomicRanges)
#' 
#' query = GRanges(rep('chr1', 5), IRanges(1:5, 1:5))
#' 
#' getPeakDistances( query )
#' 
#' @export
#' @import GenomicRanges
getPeakDistances = function( query, windowSize=10000 ){

  # expand windows
  windowRanges = query
  start(windowRanges) = start(windowRanges) - windowSize
  end(windowRanges) = end(windowRanges) + windowSize

  # get overlaps between query and expanded windows
  ovlaps = findOverlaps( windowRanges, query, ignore.strand=TRUE)
  rm(windowRanges);

  ovlapsDF = data.frame(ovlaps)
  rm(ovlaps);

  # get only top triangle because matrix is symmatrix
  ovlapsDF = ovlapsDF[with(ovlapsDF, which(queryHits >= subjectHits)),]

  # compute distances
  v1 = start(query)[ovlapsDF$queryHits] - end(query)[ovlapsDF$subjectHits]
  v2 = end(query)[ovlapsDF$queryHits] - start(query)[ovlapsDF$subjectHits]
  v3 = start(query)[ovlapsDF$queryHits] - start(query)[ovlapsDF$subjectHits]
  v4 = end(query)[ovlapsDF$queryHits] - end(query)[ovlapsDF$subjectHits]

  # for each pair, save distances
  ovlapsDF$distance = pmin(abs(v1), abs(v2), abs(v3), abs(v4))
  rm(v1); rm(v2); rm(v3); rm(v4);

  # ovlapsDF$name1 = query$names[ovlaps@from]
  # ovlapsDF$name2 = query$names[ovlaps@to]
  ovlapsDF$score1 = score(query)[ovlapsDF$queryHits]
  ovlapsDF$score2 = score(query)[ovlapsDF$subjectHits]

  attr(ovlapsDF, "N") = length(query)
  ovlapsDF
}

#' Create correlation matrix
#'
#' Create correlation matrix based on correlation between pairs of peaks
#'
#' @param query GRanges object of intervals to query
#' @param regionQuant normalized quantifications of regions in query.  Rows are features, like in limma
#' @param adjacentCount number of adjacent entries to compute correlation against
#' @param windowSize width of window in bp around each interval beyond which weight is zero
#' @param method 'adjacent': compute corr on fixed count sliding window define by adjacentCount.  "distance": compute corr for entries within windowSize bp
#' @param quiet suppress messages 
#'
#' @return for peak i and j with distance d_{i,j},
#' M[i,j] = cor( vobj$E[i,], vobj$E[j,] )
#' 
#' return sparse symmatric matrix 
#'
#' @examples
#'
#' data('decorateData')
#'
#' C = createCorrelationMatrix(simLocation, simData)
#'
#' @export
#' @importFrom Matrix sparseMatrix
createCorrelationMatrix = function( query, regionQuant, adjacentCount=500, windowSize=1e6, method = "adjacent", quiet=FALSE){

  if( nrow(regionQuant) != length(query) ){
    stop("Number of rows in regionQuant must equal entries in query")
  }

  if( method == "adjacent"){
    if( ! quiet ) cat("adjacent: ", adjacentCount, 'entries\n')
    # define pairs based on number of entries
    start(query) = seq_len(length(query))
    end(query) = seq_len(length(query))

    df_peaksMap = getPeakDistances( query, windowSize=adjacentCount )

  }else if( method == "distance"){

    if( ! quiet ) cat("distance: ", windowSize, 'bp\n')
    # define pairs based on distance
    df_peaksMap = getPeakDistances( query, windowSize=windowSize )

  }else{
    stop("method must be 'adjacent' or 'distance'")
  }

  corSubsetPairs( regionQuant, df_peaksMap$queryHits, df_peaksMap$subjectHits, silent=quiet)
}


#' Run hierarchical clustering preserving order
#' 
#' Run hierarchical clustering preserving sequential order of entries
#' 
#' @param X data matrix were *rows* are features in sequential order
#' @param gr GenomicRanges object corresponding to rows in X
#' @param alpha mixture parameter weigting correlations (alpha=0) versus chromosome distances (alpha=1)
#' 
#' @return hclust tree
#'
#' @details Use hclustgeo in ClustGeo package to generate hierarchical clustering that preserves sequential order.
#'
#' Chavent, et al. 2017. ClustGeo: an R package for hierarchical clustering with spatial constraints. arXiv:1707.03897v2 doi:10.1007/s00180-018-0791-1
#'
#' @import vegan 
#' @importFrom ClustGeo hclustgeo
#' @importFrom stats as.dist dist
#' @export
runOrderedClustering = function( X, gr, alpha=0.5 ){

  if( is.null(gr) ){
    # compute spatial restriction matrix for ordered clustering
    N = nrow(X)
    dorder = matrix(0,N,N)
    rownames(dorder) = rownames(X)
    colnames(dorder) = rownames(X)
    for(k in seq(1,N-1)){
      idx = seq(k+1, N)
      dorder[idx,k] = seq(1, (N-k))
    }
    D2 = as.dist(dorder)

  }else{
    if( nrow(X) != length(gr) ){
      stop("# rows of X must equal # of entries in gr: ", nrow(X), ' != ', length(gr))
    }

    # get squared distance between all pairs of intervals
    D2 = as.dist( outer(gr,gr, "distance") )
  }

  # compute distance matrix for data
  D1 = dist( X )
  
  treeOriginal <- hclustgeo(D1, D2, alpha=alpha)

  rm(D1)
  rm(D2)

  treeReOrder = vegan:::reorder.hclust(treeOriginal, match(treeOriginal$labels, rownames(X)))

  rm(treeOriginal)

  treeReOrder
}

setClass("epiclust", representation(clust="ANY", location="GRanges", adjacentCount="numeric", alpha="numeric", method="character", correlation="ANY"))
setClass("epiclustList", representation('ANY'))


setMethod("print", "epiclustList", function( x ){

  cat("List of hierarchical clustering for", length(x), "chromosomes:\n\n")

  counts = lapply(x, function(x) length(x@clust$labels))

  df = data.frame(chrom = names(counts), nFeatures = unlist(counts), stringsAsFactors=FALSE)
  rownames(df) = c()
  print(df, row.names=FALSE)
  cat("\n")
})

setMethod("show", "epiclustList", function( object ){
  print( object )
})

setMethod("print", "epiclust", function( x ){
  cat("Hierarchical clustering for", length(x@clust$labels), "features\n\n")
})

setMethod("show", "epiclust", function( object ){
  print( object )
})

#' Run hierarchical clustering preserving order
#' 
#' Run hierarchical clustering preserving sequential order of entries
#' 
#' @param X data matrix were *rows* are features in sequential order
#' @param gr GenomicRanges object with entries corresponding to the *rows* of X
#' @param method 'adjclust': adjacency constrind clustering.  'hclustgeo': incorporate data correlation and distance in bp
#' @param quiet suppress messages
#' @param alpha use by 'hclustgeo': mixture parameter weighing correlations (alpha=0) versus chromosome distances (alpha=1)
#' @param adjacentCount used by 'adjclust': number of adjacent entries to compute correlation against
#' 
#' @return list hclust tree, one entry for each chromosome
#'
#' @details 
#' Use adjacency constrained clustering from adjclust package:
#'
#'Alia Dehman, Christophe Ambroise and Pierre Neuvial. 2015. Performance of a blockwise approach in variable selection using linkage disequilibrium information. BMC Bioinformatics 16:148 doi.org:10.1186/s12859-015-0556-6
#'
#' Or, use hclustgeo in ClustGeo package to generate hierarchical clustering that roughly preserves sequential order.
#'
#' Chavent, et al. 2017. ClustGeo: an R package for hierarchical clustering with spatial constraints. arXiv:1707.03897v2 doi:10.1007/s00180-018-0791-1
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#' 
#' # Plot correlations and clusters in region defind by query
#' query = GRanges('chr1', IRanges(0, 1000))
#' 
#' plotDecorate( treeList, treeListClusters, query)
#'
#' @importFrom rlist list.reverse
#' @importFrom GenomicRanges seqnames
#' @importFrom adjclust adjClust
#' @importFrom methods new
#' @importFrom stats cor
#' @importFrom utils assignInNamespace
#' @export
runOrderedClusteringGenome = function( X, gr, method = c("adjclust", 'hclustgeo'), quiet=FALSE, alpha=0.5, adjacentCount=500){

  if( nrow(X) != length(gr) ){
    stop("# rows of X must equal # of entries in gr: ", nrow(X), ' != ', length(gr))
  }

  method <- match.arg(method)

  if( alpha < 0 || alpha > 1 ){
      stop("Must specify agrument: alpha must be between 0 and 1")
  }

  seqNameCounts = table(seqnames(gr))
  seqNameCounts = seqNameCounts[seqNameCounts>0]

  # chromosomes must have >- 2 entries
  tooFew = seqNameCounts[seqNameCounts<=2]
  if( length(tooFew) > 0){
    stop("Must have >=2 entires on a chromosome in order to cluster.\nThe following have too few entries: ", paste(names(tooFew), collapse=', '), "\nRemove these chromosomes before analysis")
  }

  # run last (i.e. smallest) chrom first
  chromArray = as.character(seqnames(gr)@values)

  # run analysis
  ##############

  # remove clostly sanity check in adjclust
  f = function(mat){
    NULL
  }

  assignInNamespace(x="checkCondition", value = f, ns='adjclust')

  # for each chromosome
  treePerChrom = lapply( chromArray, function(chrom){

    cat("\rEvaluating:", chrom, '          ')

    idx = which(array(GenomicRanges::seqnames(gr) == chrom))

    if(method == "adjclust"){
      C = createCorrelationMatrix( gr[idx], X[idx,], adjacentCount=adjacentCount, quiet=TRUE)
      h = min( adjacentCount, nrow(C)-1)
      fitClust = suppressMessages(adjClust( C^2, "similarity", h=h))
      # fitClust = as.hclust(fitClust)

      res = new("epiclust", clust = fitClust, location=gr[idx], adjacentCount=adjacentCount, alpha=0, method=method, correlation=C)
    }else{
      # compute tree on this chromosome
      fitClust = runOrderedClustering( X[idx,], gr[idx], alpha=alpha)

      C = cor(t(X[idx,]))

      res = new("epiclust", clust = fitClust, location=gr[idx], adjacentCount=adjacentCount, alpha=alpha, method=method, correlation=C)
    }
    rm(C)
    rm(fitClust)
    gc()
    res
  })
  cat('\n')

  names(treePerChrom) = chromArray
  
  new('epiclustList', treePerChrom )
}



#' Extract subset of data points
#' 
#' Extract subset of data points
#' 
#' @param fit epiclust or epiclustList object
#' @param query GenomicRanges object of which data points to retain
#' 
#' @return epiclust or epiclustList object
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#'
#' # extract subset of data after clustering 
#' res = getSubset( treeList, simLocation[1:10])
#' 
#' 
#' @importFrom GenomicRanges GRanges
#' @export
#' @docType methods
#' @rdname getSubset-methods
setGeneric("getSubset", function(fit, query) standardGeneric("getSubset"))

#' @importFrom GenomicRanges GRanges findOverlaps
#' @export
#' @rdname getSubset-methods
#' @aliases getSubset,epiclust,GRanges-method
setMethod("getSubset", c("epiclust", "GRanges"),
  function(fit, query) {
    # idx = NULL
    # v = match(seqnames(fit@location), seqnames(query))@values
    # if( any(!is.na(v)) > 0 ){

        # idx = which((seqnames(fit@location)@values %in% seqnames(query)@values) & (fit@location %within% query))

    # }

    ovlp = suppressWarnings(findOverlaps(fit@location, query))
    idx = ovlp@from

    if( length(idx) == 0){
      fitNew = NULL
    }else{
      fitNew = fit

      N = length(idx)

      if( N == 0 ){
        stop("Cannot select empty subset of entries")
      }
      idx2 = idx[seq_len(length(idx)-1)]
      fitNew@clust$merge = fitNew@clust$merge[idx2,]
      fitNew@clust$height = fitNew@clust$height[idx2]
      fitNew@clust$order = fitNew@clust$order[idx]
      fitNew@clust$labels = fitNew@clust$labels[idx]

      fitNew@location = fitNew@location[idx]

      fitNew@correlation = fitNew@correlation[idx,idx]
    }
    fitNew
  }
)

#' @importFrom GenomicRanges GRanges
#' @export
#' @rdname getSubset-methods
#' @aliases getSubset,epiclustList,GRanges-method
setMethod("getSubset", c("epiclustList", "GRanges"),
   function(fit, query){
    res = lapply(fit, function(x) getSubset(x, query))
    Filter(Negate(is.null), res)
  }
)



#' Create cluster from list of hclust objects
#'
#' Create cluster from list of hclust objects
#'
#' @param treeList list of hclust objects 
#' @param method 'capushe': slope heuristic. 'bstick': broken stick. 'meanClusterSize': create clusters based on target mean value. 
#' @param meanClusterSize select target mean cluster size
#' @param pct minimum percentage of points for the plateau selection in capushe selection.
#'
#' @return  Convert hierarchical clustering into discrete clusters based on selection criteria method
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#' 
#' # Plot correlations and clusters in region defind by query
#' query = GRanges('chr1', IRanges(0, 1000))
#' 
#' plotDecorate( treeList, treeListClusters, query)
#'
#' @import adjclust
#' @importFrom stats cutree
#' @export
createClusters = function(treeList, method = c("capushe", "bstick", "meanClusterSize"), meanClusterSize=50, pct=0.15) {
  
  method <- match.arg(method)

  cat("Method:", method, '\n')

  if( method == "meanClusterSize"){
    # get counts on each chromosome
    n_features = vapply(treeList, function(x){
       length(x@clust$order)
    }, FUN.VALUE=numeric(1))
    # combine counts across chromosomes
    n_features_total = sum(n_features)

    n_total_clusters = round(n_features_total / meanClusterSize)

    # cut tree to get average of meanClusterSize
    res = lapply(treeList, function(x){
      N = length(x@clust$order)
      frac = N / n_features_total
      n_clust = round( frac*n_total_clusters )
      cutree(x@clust, k=n_clust)
    })
    names(res) = names(treeList)

  }else{

    if( treeList[[1]]@method == 'hclustgeo' ){
      stop("Method '", method, "' is only compatible with result of\nrunOrderedClusteringGenome(...,method ='adjclust')")
    }

    res = lapply(treeList, function(x){
      select( x@clust, type=method, pct=pct)
    })
    names(res) = names(treeList)     
  }
  new('epiclustDiscreteList',res)
}

# setClass("epiclustDiscrete", representation("array"))

# setMethod("print", "epiclustDiscrete", function( x ){
#   nFeatures = length(x)
#   nClusters = length(unique(x))
#   cat(nFeatures, "features in", nClusters, "discrete clusters\n\n" )
# })

# setMethod("show", "epiclustDiscrete", function( object ){
#   print( object )
# })


setClass("epiclustDiscreteList", representation('list'))

setMethod("print", "epiclustDiscreteList", function( x ){
  for( chrom in names(x) ){
    nFeatures = length(x[[chrom]])
    nClusters = length(unique(x[[chrom]]))
    cat(paste0(chrom, ':'),nFeatures, "features in", nClusters, "discrete clusters\n\n" )
  }
})

setMethod("show", "epiclustDiscreteList", function( object ){
  print( object )
})






#' Simulated data to show correlation clustering
#'
#' @docType data
#' @keywords decorateData
#' @name decorateData
#' @usage data(decorateData)
"simData"

#' Simulated data to show correlation clustering. GRanges object indicating genomic position of data in rows of simData
#'
#' @docType data
#' @keywords decorateData
#' @name decorateData
#' @usage data(decorateData)
"simLocation"




#' Count clusters on each chromosome
#'
#' Count clusters on each chromosome
#'
#' @param treeListClusters from createClusters()
#'
#' @return  count number of clusters on each chromsome
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#' 
#' # Count clusters on each chromsome
#' countClusters( treeListClusters )
#'
#' @export
countClusters = function(treeListClusters){
 vapply(treeListClusters, function(x) length(unique(x)), numeric(1))
}


#' Find which cluster a peak is in
#'
#' Find which cluster a peak is in
#'
#' @param treeListClusters from createClusters()
#' @param id name of query peak
#'
#' @return data.frame of chromosome and cluster
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#' 
#' # Find chromsome and cluster of peak_204
#' whichCluster( treeListClusters, 'peak_204')
#'
#' @export
whichCluster = function(treeListClusters, id){
  # find cluster based on anem
  res = lapply(names(treeListClusters), function(chrom){
    x = treeListClusters[[chrom]]
    idx = match(id, names(x))
    data.frame( chrom, id, cluster=x[idx], stringsAsFactors=FALSE)
    })
  res = do.call('rbind', res)

  res[!is.na(res$cluster),]
}


# whichCluster = function(treeListClusters, id){
#   # find cluster based on anem
#   res = vapply(treeListClusters, function(x){
#     idx = (names(x) == id)
#     if( sum(idx) == 0){
#       res = NA
#     }else{
#       res = x[idx]
#     }
#     res}, numeric(1))

#   res = res[!is.na(res)]

#   if( length(res) == 0){
#     chrom = cluster = NA
#   }else{
#     chrom = names(res)
#     cluster = res[1]
#   }
#   res = data.frame(chrom=chrom, cluster=cluster, stringsAsFactors=FALSE)
#   rownames(res) = c()
#   res
# }


# whichCluster(treeListClusters, 'Peak_143742')
# whichCluster2(treeListClusters, 'Peak_143742')

#' Get feature names in selected cluster
#'
#' Get feature names in selected cluster given chrom and clusterid 
#'
#' @param treeListClusters from createClusters()
#' @param chrom chromosome name of cluster
#' @param clustID cluster identifier
#'
#' @return array of feature names
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#' 
#' # Find chromsome and cluster of peak_204
#' getFeaturesInCluster( treeListClusters, "chr1", 3)
#'
#' @export
getFeaturesInCluster = function( treeListClusters, chrom, clustID){

  idx = which( treeListClusters[[chrom]] == clustID )

  clsts = treeListClusters[[chrom]][idx]

  names(clsts)
}



#' Compute scores for each cluster
#'
#' For each cluster compute summary statistics for the cluster to measure how strong the correlation structure is.  Clusters with weak correlation structure can be dropped from downstream analysis.
#'
#' @param treeList list of hclust objects 
#' @param treeListClusters from createClusters()
#' @param BPPARAM parameters for parallel evaluation
#' 
#' @details For each cluster, extract the correlation matrix and return the mean absolute correlation; the 75th, 90th and 95th quantile absolute correlation, and LEF, the leading eigen-value fraction which is the fraction of variance explained by the leading eigen value of the matrix abs(C).
#'
#' @return for all pairs of peaks within windowSize, report distance 
#' 
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#'
#' # Evaluate score for each cluster
#' clstScore = scoreClusters(treeList, treeListClusters )
#' 
#' @export
#' @importFrom stats quantile
#' @importFrom progress progress_bar
scoreClusters = function(treeList, treeListClusters, BPPARAM=SerialParam()){

  # get data.frame of cluster names per chromosomes
  clCount = countClusters( treeListClusters )

  df_count = lapply( seq_len(length(clCount)), function(i){
    x = clCount[i]
    data.frame(chrom=rep(names(x), x), cluster=seq(x), stringsAsFactors=FALSE)
    })
  df_count = do.call("rbind", df_count)

  # evaluate for each cluster
  .score_clust = function(i, treeList, treeListClusters, df_count, pb){
    
    pb$update( i / nrow(df_count) )      

    # get chrom and clsuter ID for this clsuter
    chrom = df_count$chrom[i] 
    clustID = df_count$cluster[i] 

    cl = treeListClusters[[chrom]] 

    peakIDs = names(cl[cl==clustID])

    N = length(peakIDs)

    if( N > 1){
      # Get correlation matrix
      C = treeList[[chrom]]@correlation[peakIDs,peakIDs]

      # eigen analysis
      dcmp = eigen(abs(C), symmetric=TRUE, only.values = TRUE)
      LEF = dcmp$values[1] / sum(dcmp$values)

      # get correlation values
      diag(C) = -Inf
      cor_values = C@x[is.finite(C@x)]
    }else{
      cor_values = 1
      LEF = 1
    }
    # compute statistics
    data.frame( chrom = chrom, 
        cluster = clustID, N=N,
        mean_abs_corr = mean(abs(cor_values)), 
        quantile75 = quantile(abs(cor_values), .75), 
        quantile90 = quantile(abs(cor_values), .90), 
        quantile95 = quantile(abs(cor_values), .95),
        LEF=LEF, stringsAsFactors=FALSE)
  }
  
  # .score_clust(2477, treeList, treeListClusters, df_count, pb)
  
  pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta", total = sum(clCount), width= 60, clear=FALSE)

  cat("Evaluating strength of each cluster\n\n")

  # Evaluate score for each cluster
  clustScores = lapply( seq_len(nrow(df_count)), .score_clust, treeList, treeListClusters, df_count, pb)

  # format for return
  clustScores = do.call("rbind", clustScores)
  rownames(clustScores) = c()
  clustScores
}

#' Extract subset of clusters
#'
#' Extract subset of clusters based on entries in chroms and clusters
#'
#' @param treeListClusters from createClusters()
#' @param chroms array of chromosome names for the clusters to retain
#' @param clusters array of cluster names for the clusters to retain
#' 
#' @details Each cluster has a unique identifier: c(chroms[i], clusters[i]).  Return epiclustDiscreteList object like treeListClusters but retaining on the specified clusters.
#'
#' @return epiclustDiscreteList of specified clusters
#' 
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clusters
#' treeListClusters = createClusters( treeList )
#'
#' # Evaluate score for each cluster
#' clstScore = scoreClusters(treeList, treeListClusters )
#' 
#' # Retain clusters that pass this criteria
#' clustInclude = clstScore[clstScore$LEF > 0.2,]
#' 
#' # get retained clusters
#' treeListClusters_filter = filterClusters( treeListClusters, clustInclude$chrom, clustInclude$cluster)
#' 
#' @export
filterClusters = function( treeListClusters, chroms, clusters){

  if( ! is(treeListClusters, 'epiclustDiscreteList') ){
    stop("treeListClusters must be of class 'epiclustDiscreteList'")
  }

  if( length(chroms) < 1 ){
    stop("length of chroms must be >= 1: ", length(chroms), ' < 1')
  }
  if( length(chroms) != length(clusters) ){
    stop("chroms and clusters must have same length: ",  length(chroms), ' != ', length(clusters))
  }

  # convert to data.frame
  df = data.frame(chrom = chroms, cluster=clusters, stringsAsFactors=FALSE)

  # extract from treeListClusters entries corresponding to chroms and clusters
  resList = lapply( unique(df$chrom), function(chrom){
    # get indecies in df
    idx = which(chrom == df$chrom)

    # get indeces in treeListClusters
    idx2 = treeListClusters[[chrom]] %in% df$cluster[idx]

    # get clusters in treeListClusters
    treeListClusters[[chrom]][idx2]
    })
  names(resList) = unique(df$chrom)

  new('epiclustDiscreteList', resList)
}


#' Evaluate the decay of correlation versus distance between features
#'
#' For pairs of features evaluate the physical distance and the correlation
#'
#' @param treeList list of hclust objects
#' @param gr GenomicRanges object corresponding to features clustered in treeList
#'
#' @return a data.frame of distance and correlation value for all pairs of features already evalauted in treeList.  Note that runOrderedClusteringGenome() that returns treeList only evalutes correlation between a specified number of adjacent peaks 
#' 
#' @examples
#' library(GenomicRanges)
#' library(ggplot2)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clustering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#'
#' # Evaluate how correlation between features decays with distance
#' dfDist = evaluateCorrDecay( treeList, simLocation )
#' 
#' # make plot
#' ggplot(dfDist, aes(distance, abs(correlation))) + theme_bw(17) + theme(aspect.ratio=1) + geom_smooth(se=FALSE)
#'
#' @import methods
#' @import Matrix
#' @import GenomicRanges
#' @importFrom data.table data.table
#' @export
evaluateCorrDecay = function( treeList, gr){

  distList = lapply( levels(seqnames(gr)), function( chrom ){
    # get GenomicRange for this chromosome
    gRange = gr[seqnames(gr) == chrom]

    # get correlation matrix for this chromosome
    C = treeList[[chrom]]@correlation

    # convert format of sparse correlation matrix
    A = as(C, 'TsparseMatrix')

    # extract index and correlation value for each non-zero pair
    dfDist = data.frame( chrom=chrom,
      feature_i=rownames(A)[A@i+1], i=A@i+1, 
                         feature_j=rownames(A)[A@j+1], j=A@j+1, 
                         correlation = A@x,
                         stringsAsFactors=FALSE ) 

    # compute chromsomal distance between pairs with non-zero correlation
    dfDist$distance = distance(gRange[dfDist$i],gRange[dfDist$j])
    dfDist
  })

  feature_i = feature_j = NA
  
  dfDist = do.call("rbind", distList)
  dfDist = data.table(dfDist)
  dfDist = dfDist[feature_i!=feature_j,]

  dfDist[,c("chrom", "feature_i", "feature_j", "correlation", "distance")]
}


