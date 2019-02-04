
#' Plot triangle of correlation matrix
#' 
#' Plot lower triangle of correlation matrix
#' 
#' @param C correlation matrix
#' @param size plotting argument to geom_point()
#' @param stroke plotting argument to geom_point()
#' @param cols array of two colors for gradient
#' 
#' @return ggplot2 plot of correlation matrix
#'
#' @details Adjust size and stroke of points in the plot to fix look of plot depending on dimensions
#'
#' @examples
#' N = 1000
#' p = 100
#' X = matrix(rnorm(N*p), N,p)
#' C = cor(X)
#' plotCorrTriangle( C )
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom methods is
#' @export 
plotCorrTriangle = function(C, size=1, stroke=1.5, cols=c("blue", "white","red")){

  # to pass R CMD check
  x = y = value = NULL

  if( is(C, "matrix")){
    C[upper.tri(C, diag=TRUE)] = NA
    df = melt(C)
    df = df[!is.na(df$value),]    
    colnames(df)[1] = "x"
    colnames(df)[2] = "y"
    df$x = as.numeric(df$x)
    df$y = as.numeric(df$y)
  }else{
    # convert sparse matrix
    df = summary(C)
    colnames(df) = c("y", "x", "value")
    df = df[df$x!=df$y,]
    df = rbind(df, c(1, nrow(C), NA))
  }

  rotate <- function(df, degree) {
    dfr <- df
    degree <- pi * degree / 180
    l <- sqrt(df$x^2 + df$y^2)
    teta <- atan(df$y / df$x)
    dfr$x <- (l * cos(teta - degree))
    dfr$y <- (l * sin(teta - degree))
    return(dfr)
  }

  df_rotate = rotate( df, degree=45 )

  # set aspect ratio to be sqrt(1/2) so that each box is square
  # due to Pythagorean ratio
  fig = ggplot(df_rotate, aes(x, y)) + geom_point(aes(color=value, fill=value), size=size, stroke=stroke, shape=18) + scale_color_gradientn(name = "Correlation", colours=cols, limits=c(-1,1), na.value="white") + scale_fill_gradientn(name = "Correlation", colours=cols, limits=c(-1,1), na.value="white")  + theme_void() + theme(aspect.ratio=1/sqrt(2), plot.title = element_text(hjust = 0.5), legend.position="bottom") #+ xlim(range(df_rotate$x)) # 

  # get range along x-axis
  attr(fig, 'xrange') = range(df_rotate$x)
  fig
}


#' Plot cluster segments
#'
#' Plot bar of segments showing clusters
#'
#' @param clusterValues array of names of cluster for each entry
#'
#' @examples
#' plotClusterSegments(c(rep(1, 5), rep(2,2), rep(3, 4)))
#'
#' @return ggplot2 of cluster assignments
#' @import ggplot2
#' @export
plotClusterSegments = function( clusterValues ){

  dfclust = data.frame(cluster=as.character(clusterValues), stringsAsFactors=FALSE)
  dfclust$x = 1
  dfclust$y = seq_len(nrow(dfclust)) *sqrt(2)

  # to pass R CMD check
  x = y = cluster = NULL

  ggplot(dfclust, aes(x, y, fill=cluster)) + geom_tile()  + coord_flip() + theme_void() + theme(legend.position="none", aspect.ratio=0.05)
}



#' Plot decorate analysis
#'
#' Plot decorate analysis for clusters and correlations
#'
#' @param treeList hierarchical clustering of each chromosome from runOrderedClusteringGenome()
#' @param treeListClusters assign regions to clusters after cutting tree with createClusters()
#' @param query GRanges object indiecating region to plot
#' @param size plotting argument to geom_point() in correlation triangle
#' @param stroke plotting argument to geom_point() in correlation triangle
#' @param cols array of two colors for gradient in correlation triangle
#' @param plotTree show tree from hierachical clustering
#'
#' @return ggplot2 of cluster assignments and correlation between peaks
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clsutering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clutsers
#' treeListClusters = createClusters( treeList )
#' 
#' # Plot correlations and clusters in region defind by query
#' query = GRanges('chr1', IRanges(0, 1000))
#' 
#' plotDecorate( treeList, treeListClusters, query)
#'
#' @import ggplot2
#' @import grid
#' @import Biobase
#' @importFrom labeling extended
#' @importFrom stats as.hclust
#' @importFrom ggdendro ggdendrogram
#' @importFrom adjclust correct
# @import BiocGenerics
#' @export
plotDecorate = function( treeList, treeListClusters, query, size=1, stroke=1.5, cols=c("white","red"), plotTree=TRUE){

  if( length(query) > 1){
    stop("Can only query one interval")
  }

  fit = getSubset( treeList, query )

  if( length(fit) == 0){
    stop("No data left after filtering by query region")
  }

  chrom = names(fit)

  clst = treeListClusters[[chrom]]
  idx = names(clst) %in% fit[[1]]@clust$labels
  clst = clst[idx]

  if( any(!idx) & plotTree){
    cat("Cannot plot tree when features are dropped\n")
    plotTree = FALSE
  }

  # plot correlation matix
  fig1 = plotCorrTriangle( fit[[1]]@correlation, size=size, stroke=stroke, cols=cols )

  fig2 = plotClusterSegments( clst )

  # scale
  rng = range(fit[[1]]@location)
  v1 = start(rng) 
  v2 = end(rng) 

  m = 5
  breaks = extended(v1, v2, m)

  # df = data.frame(x=c(v1,v2), y=c(0,0))
  df = data.frame(x=c(min(breaks[1], v1),max(v2, breaks[m])), y=c(0,0))

  # to pass R CMD check
  x = y = NULL

  fig_Scale = ggplot(df, aes(x,y)) + theme_bw() + 
    theme(aspect.ratio=0.1,
      plot.title = element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y=element_blank()) + geom_line() + xlab('') + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + xlab('') + scale_x_continuous(limits=range(df$x), breaks=breaks) 

  if( plotTree ){

    treeTmp = fit[[1]]@clust 

    # if tree is chac from adjclust, use correct() to get positive branch lengths
    if( is(treeTmp, 'chac') ){
      treeTmp = as.hclust(correct( treeTmp ))
    }

    figTree = ggdendrogram(treeTmp, labels=FALSE, leaf_labels=FALSE) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())

    figCombine = rbind(ggplotGrob(figTree), ggplotGrob(fig_Scale), ggplotGrob(fig2), ggplotGrob(fig1), size="last" )
  }else{
    figCombine = rbind( ggplotGrob(fig_Scale), ggplotGrob(fig2), ggplotGrob(fig1), size="last")
  }
   
  grid.newpage()
  grid.draw( figCombine )
}



#' Plot two correlation matrices together
#'
#' Combined plot of correlation matricies from cases and controls
#'
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param peakIDs feature names to extract from rows of epiSignal
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param size size of text
#' @param cols array of 3 color values
#'
#' @return ggplot2 of combined correlation matrix
#'
#' @examples
#' library(GenomicRanges)
#' 
#' data('decorateData')
#' 
#' # Evaluate hierarchical clsutering
#' treeList = runOrderedClusteringGenome( simData, simLocation ) 
#' 
#' # Choose cutoffs and return clutsers
#' treeListClusters = createClusters( treeList )
#'
#' # Simulate variable to split dataset by
#' set.seed(1)
#' metadata = data.frame( Disease = factor(sample(0:1, ncol(simData), replace=TRUE)))
#'
#' # get peak ID's from chr1, cluster 1
#' peakIDs = getFeaturesInCluster( treeListClusters, "chr1", 1)
#'
#' # plot comparison of correlation matrices for peaks in peakIDs
#' #  where data is subset by metadata$Disease
#' plotCompareCorr( simData, peakIDs, metadata$Disease) + ggtitle("chr1: cluster 1")
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotCompareCorr = function(epiSignal, peakIDs, testVariable, size=5, cols=c("blue", "white","red")){

  if( length(testVariable) != ncol(epiSignal) ){
    stop("Number of columns in epiSignal must equal number of entries in testVariable")
  }
  if( ! is(testVariable, 'factor') || nlevels(testVariable) != 2 ){
    stop("Entries in testVariable must be a factor with two levels")
  }

  # Divide data into two sets
  set1 = which(testVariable == levels(testVariable)[1])
  set2 = which(testVariable == levels(testVariable)[2])

  # get features
  Var1 = Var2 = value = NA
  Y1 = scale(t(epiSignal[peakIDs, set1]))
  Y2 = scale(t(epiSignal[peakIDs, set2]))

  # evaluate correlation
  C1 = cor(Y1)
  C2 = cor(Y2)

  # create combined correlation matrix
  C = C1
  C[lower.tri(C)] = C2[lower.tri(C2)]
  diag(C) = NA

  # convert matrix to data.frame to plot
  df = melt( C )
  N = nrow( C )

  ggplot(df, aes(Var1, Var2)) + geom_tile(aes(color=value, fill=value)) + scale_color_gradientn(name = "Correlation", colours=cols, limits=c(-1,1), na.value="grey") + scale_fill_gradientn(name = "Correlation", colours=cols, limits=c(-1,1), na.value="grey")  + theme_void() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="bottom") + annotate(geom="text", x=c(0.1*N,0.9*N), y=c(0.9*N,0.1*N), label=levels(testVariable),size=size)
}







