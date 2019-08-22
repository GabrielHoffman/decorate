
#' Evaluate the decay of correlation versus distance between features
#'
#' For pairs of features evaluate the physical distance and the correlation
#'
#' @param treeList list of hclust objects
#' @param gr GenomicRanges object corresponding to features clustered in treeList
#' @param chromArray Use this only this set of chromosmes.  Can substantially reduce memory usage
#' @param verbose show progress
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
#' plotCorrDecay( dfDist )
#'
#' @import methods
#' @import Matrix
#' @import GenomicRanges
#' @importFrom data.table data.table
#' @export
evaluateCorrDecay = function( treeList, gr, chromArray=seqlevels(gr), verbose=TRUE){
  
  # Drop empty chromsomes
  chrNameCount = table(seqnames(gr))
  gr = dropSeqlevels( gr, names(chrNameCount[chrNameCount==0]))

  # chromArray should only contain chromsome with at least 1 peak
  chromArray = chromArray[chromArray%in% seqlevels(gr)]

  distList = list()
  for( chrom in seqlevels(gr)[seqlevels(gr) %in% chromArray] ){

    if( verbose ){
      cat(paste0("\r", chrom, '      '))
    }

    # get GenomicRange for this chromosome
    gRange = gr[seqnames(gr) == chrom]

    # get correlation matrix for this chromosome
    # convert format of sparse correlation matrix
    A = as(treeList[[chrom]]@correlation, 'TsparseMatrix')

    # extract index and correlation value for each non-zero pair
    dfDist = data.frame( chrom=chrom,
      feature_i=rownames(A)[A@i+1], i=A@i+1, 
                         feature_j=rownames(A)[A@j+1], j=A@j+1, 
                         correlation = A@x,
                         stringsAsFactors=FALSE ) 
    rm(A)

    # compute chromsomal distance between pairs with non-zero correlation
    dfDist$distance = distance(gRange[dfDist$i],gRange[dfDist$j])
    distList[[chrom]] = dfDist
  }

  if( verbose ){
    cat("\n")
  }

  feature_i = feature_j = NA
  
  dfDist = do.call("rbind", distList)
  rm(distList)
  dfDist = data.table(dfDist)
  dfDist = dfDist[feature_i!=feature_j,]

  dfDist[,c("chrom", "feature_i", "feature_j", "correlation", "distance")]
}

#' Plot correlation delay
#'
#' Plot correlation delay using subsampling
#'
#' @param dfDist data.frame of distance and correlation from from evaluateCorrDecay() 
#' @param n number of points to sample per bin
#' @param nbins number of bins on the x-axis
#' @param method on show either R or Rsq on y-axis
#' @param xlim min and max values for x-axis
#'
#' @details Plot correlation versus log10 distance.  Sample equal number of points for each bin along the x-axis.
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
#' plotCorrDecay( dfDist )
#'
#' @import ggplot2
#' @export
plotCorrDecay = function( dfDist, n = 20000, nbins=1000, method = c("R", "Rsq"), xlim=c(10,1e6) ){  

  method = match.arg(method)

  if( method == 'R' ){
    fig = ggplot_by_sampling( log10(dfDist$distance), dfDist$correlation, n, nbins ) + ylim(-1,1) + ylab("Correlation") 
  }else{
     fig = ggplot_by_sampling( log10(dfDist$distance), dfDist$correlation^2, n, nbins ) + ylim(0,1) + ylab("Squared Correlation") 
  }

  breaks = c(0, 1, 2,3,4, 5, 6, 7)
  labels = c('1bp', '10 bp', '100bp', '1 kb', '10 kb', '100 kb', '1 Mb', '10 Mb')

  idx = (10^breaks >= xlim[1]) & (10^breaks <= xlim[2])

  fig + geom_point(size=.1) + theme_bw(17) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + geom_smooth() + ggtitle("Correlation versus distance between features") + xlab(bquote(Distance~between~features~(log[10]~scale))) + scale_x_continuous(breaks=breaks[idx], labels=labels[idx], limits=log10(xlim))
}


#' Plot by subsampling in each bin 
#'
#' Plot by subsampling in each bin in x-axis
#'
#' @param x x values
#' @param y y values
#' @param N number of samples
#' @param nbins number of bins on the x-axis
#'
ggplot_by_sampling = function(x, y, N, nbins=1000){

  # define tp pass R CMD check
  cluster = prob = NA

  # convert to data.table
  df = data.table(x,y)

  # omit rows with NA
  df = df[complete.cases(df*0)]

  # order by distance
  df = df[order(x),]

  # create bins
  bins = seq(min(df$x), max(df$x), length.out=nbins)

  # assign each point to a cluster based on distance
  df[,cluster:=as.character(findInterval(x, bins))]

  # count number of entries in each cluster
  tab = df[,table(cluster)]

  # merge to assign counts for each cluster
  df = merge(df, data.table(tab), by="cluster")

  # set probability proportion to 1/N
  df[,prob:=1/as.numeric(N)]

  # sample indeces
  idx = sample.int(nrow(df), N, replace=TRUE, prob = df$prob )

  # make plot using only unique indieces
  ggplot(df[unique(idx),], aes(x, y))
}
