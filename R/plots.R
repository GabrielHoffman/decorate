
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
plotCorrTriangle = function(C, size=1, stroke=1.5, cols=c("white","red")){

  # to pass R CMD check
  x = y = value = NULL

  if( is(C) == "matrix"){
    C[upper.tri(C, diag=TRUE)] = NA
    df = melt(abs(C))
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
  fig = ggplot(df_rotate, aes(x, y)) + geom_point(aes(color=value, fill=value), size=size, stroke=stroke, shape=18) + scale_color_gradientn(name = "Correlation", colours=cols, limits=c(0,1), na.value="white") + scale_fill_gradientn(name = "Correlation", colours=cols, limits=c(0,1), na.value="white")  + theme_void() + theme(aspect.ratio=1/sqrt(2), plot.title = element_text(hjust = 0.5), legend.position="bottom") #+ xlim(range(df_rotate$x)) # 

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
  dfclust$y = 1:nrow(dfclust) *sqrt(2)

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
#'
#' @return ggplot2 of cluster assignments and correlation between peaks
#'
#' @examples
#' library(decorate)
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
#' @export
plotDecorate = function( treeList, treeListClusters, query){

  if( length(query) > 1){
    stop("Can only query one interval")
  }

  fit = getSubset( treeList, query )

  if( length(fit) == 0){
    stop("No data left after filtering by query region")
  }

  chrom = names(fit)

  clst = treeListClusters[[chrom]]
  clst = clst[names(clst) %in% fit[[1]]@clust$labels]

  # plot correlation matix
  fig1 = plotCorrTriangle( fit[[1]]@correlation, cols=c("white", "red") )

  fig2 = plotClusterSegments( clst )

  # scale
  rng = range(fit[[1]]@location)
  v1 = start(rng) 
  v2 = end(rng) 

  df = data.frame(x=c(v1,v2), y=c(0,0))

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
      axis.ticks.y=element_blank()) + geom_line() + scale_x_continuous(limits=c(v1,v2)) + xlab('') + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + xlab('') 

  figCombine = rbind(ggplotGrob(fig_Scale), ggplotGrob(fig2), ggplotGrob(fig1), size="last" )
   
  grid.newpage()
  grid.draw( figCombine )
}





