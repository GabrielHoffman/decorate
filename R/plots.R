
#' Plot triangle of correlation matrix
#' 
#' Plot lower triangle of correlation matrix
#' 
#' @param C correlation matrix
#' @param size plotting argument to geom_point()
#' @param stroke plotting argument to geom_point()
#' @param cols array of two colors for gradient
#' @param absCorr show absolute correlations
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
plotCorrTriangle = function(C, size=1, stroke=1.5, cols=c("blue", "white","red"), absCorr=FALSE){

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

  if( absCorr ){
    limits=c(0,1)
    cols = cols[2:3]
    df$value = abs( df$value )
  }else{
    limits=c(-1,1)  
  }

  df_rotate = rotate( df, degree=45 )

  # set aspect ratio to be sqrt(1/2) so that each box is square
  # due to Pythagorean ratio
  fig = ggplot(df_rotate, aes(x, y)) + geom_point(aes(color=value, fill=value), size=size, stroke=stroke, shape=18) + scale_color_gradientn(name = "Correlation", colours=cols, limits=limits, na.value="white") + scale_fill_gradientn(name = "Correlation", colours=cols, limits=limits, na.value="white")  + theme_void() + theme(aspect.ratio=1/sqrt(2), plot.title = element_text(hjust = 0.5), legend.position="bottom") #+ xlim(range(df_rotate$x)) # 

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

  ggplot(dfclust, aes(x, y, fill=cluster)) + geom_tile()  + coord_flip() + theme_void() + theme(legend.position="none", aspect.ratio=0.05) + scale_fill_discrete(na.value = 'white')
}


# ' Plot decorate analysis
# '
# ' Plot decorate analysis for clusters and correlations
# '
# ' @param treeList hierarchical clustering of each chromosome from runOrderedClusteringGenome()
# ' @param treeListClusters assign regions to clusters after cutting tree with createClusters()
# ' @param query GRanges object indiecating region to plot
# ' @param size plotting argument to geom_point() in correlation triangle
# ' @param stroke plotting argument to geom_point() in correlation triangle
# ' @param cols array of two colors for gradient in correlation triangle
# ' @param plotTree show tree from hierachical clustering
# ' @param absCorr show absolute correlations
# '
# ' @return ggplot2 of cluster assignments and correlation between peaks
# '
# ' @examples
# ' library(GenomicRanges)
# ' 
# ' data('decorateData')
# ' 
# ' # Evaluate hierarchical clsutering
# ' treeList = runOrderedClusteringGenome( simData, simLocation ) 
# ' 
# ' # Choose cutoffs and return clusters
# ' treeListClusters = createClusters( treeList )
# ' 
# ' # Plot correlations and clusters in region defind by query
# ' query = GRanges('chr1', IRanges(0, 1000))
# ' 
# ' plotDecorate( treeList, treeListClusters, query)
# '
# ' @import ggplot2
# ' @import grid
# ' @import Biobase
# ' @importFrom labeling extended
# ' @importFrom stats as.hclust
# ' @importFrom ggdendro ggdendrogram
# ' @importFrom adjclust correct
# @import BiocGenerics
# ' @export
# plotDecorate = function( treeList, treeListClusters, query, size=1, stroke=1.5, cols=c("blue", "white","red"), plotTree=TRUE, absCorr = TRUE){

#   if( length(query) > 1){
#     stop("Can only query one interval")
#   }

#   fit = getSubset( treeList, query )

#   if( length(fit) == 0){
#     stop("No data left after filtering by query region")
#   }

#   chrom = names(fit)

#   clst = treeListClusters[[chrom]]
#   idx = names(clst) %in% fit[[1]]@clust$labels
#   clst = clst[idx]

#   # create new set of clusters including non-clusters
#   clstComplete = rep(NA,nrow(fit[[1]]@correlation))
#   names(clstComplete) = rownames(fit[[1]]@correlation)

#   idx2 = match(names(clst), names(clstComplete))
#   clstComplete[idx2] = clst

#   if( any(!idx) & plotTree){
#     cat("Cannot plot tree when features are dropped\n")
#     plotTree = FALSE
#   }

#   # plot correlation matix
#   fig1 = plotCorrTriangle( fit[[1]]@correlation, size=size, stroke=stroke, cols=cols, absCorr=absCorr )


#   fig1 = plotCorrTriangle2( fit[[1]]@correlation )

#   fig2 = plotClusterSegments( clstComplete )

#   # scale
#   # rng = range(fit[[1]]@location)
#   # v1 = start(rng) 
#   # v2 = end(rng) 

#   # m = 5
#   # breaks = extended(v1, v2, m)

#   # # df = data.frame(x=c(v1,v2), y=c(0,0))
#   # df = data.frame(x=c(min(breaks[1], v1),max(v2, breaks[m])), y=c(0,0))

#   # # to pass R CMD check
#   # x = y = NULL

#   # fig_Scale = ggplot(df, aes(x,y)) + theme_bw() + 
#   #   theme(aspect.ratio=0.1,
#   #     plot.title = element_blank(),
#   #     axis.title.y=element_blank(),
#   #     axis.text.y=element_blank(), 
#   #     panel.grid.major = element_blank(), 
#   #     panel.grid.minor = element_blank(),
#   #     panel.border = element_blank(),
#   #     axis.ticks.y=element_blank()) + geom_line() + xlab('') + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + xlab('') + scale_x_continuous(limits=range(df$x), breaks=breaks) 

#   if( plotTree ){

#     treeTmp = fit[[1]]@clust 

#     # if tree is chac from adjclust, use correct() to get positive branch lengths
#     if( is(treeTmp, 'chac') ){
#       treeTmp = as.hclust(correct( treeTmp ))
#     }

#     figTree = ggdendrogram(treeTmp, labels=FALSE, leaf_labels=FALSE) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())

#     figCombine = rbind(ggplotGrob(figTree), ggplotGrob(fig2), ggplotGrob(fig1), size="last" )
#     #, ggplotGrob(fig_Scale)
#   }else{
#     figCombine = rbind( ggplotGrob(fig2), ggplotGrob(fig1$LDheatmapGrob), size="last")
#     #, ggplotGrob(fig_Scale)
#   }
   
#   grid.newpage()
#   grid.draw( figCombine )
# }





#' Convert correlation matrix into triangle plot
#'
#' Adapted from LDheatmap
#'
#' @param nrow nrow(C)
#' @param ncol ncol(C)
#' @param entries of C convertned to color
#' @param name name of the plot
#' @param byrow process C by row
#' @param cutoff retain correlations up to cutoff features away 
#'
makeImageRect <- function(nrow, ncol, cols, name, byrow=TRUE, cutoff=100) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  # Creates coordinate pairs, repeating either column numbers (if byrow = TRUE) or row numbers (if byrow = FALSE) to force that type of fill
  if(byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }

  i = (right-top)^2 < cutoff/nrow

  rectGrob(x=right[i], y=top[i], 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols[i]),
           name=name)
}


addLegend <- function(color, vp){
  ImageRect<- makeImageRect(2,length(color), cols=c(rep(NA,length(color)),color[length(color):1]), "colorKey", cutoff=1000)

  keyVP <- viewport(x=.94, y=.03, height=.05, width=.2, just=c("right","bottom"), name="keyVP")
  
  ttt<-expression(paste(R^2))  
  title<-textGrob(ttt, x=0.5, y=1.25, name="title", gp=gpar(cex=0.8))
  
  #Adding labels to the color key
  labels<-textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="labels")
  
  #Drawing ticks at the bottom axis of the color key
  ticks<-segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="ticks")
  
  #Drawing a box around the color key
  box <- linesGrob(x=c(0,0,1,1,0), y=c(0.5,1,1,0.5,0.5), name="box")
  
  key <- gTree(children=gList(ImageRect, title, labels, ticks, box), name = "Key", vp=keyVP)
  key
}

#' Plot decorate analysis
#' 
#' Plot decorate analysis for clusters and correlations
#'
#' @param treeList hierarchical clustering of each chromosome from runOrderedClusteringGenome()
#' @param treeListClusters assign regions to clusters after cutting tree with createClusters()
#' @param query GRanges object indiecating region to plot
#' @param plotTree show tree from hierachical clustering
#' @param windowWidth for each feature plot correlation with windowWidth adjacent features 
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
#' # Choose cutoffs and return clusters
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
#' @import GenomicRanges
#' @importFrom labeling extended
#' @importFrom stats as.hclust
#' @importFrom ggdendro ggdendrogram
#' @importFrom reshape2 melt
#' @importFrom adjclust correct
#' @export
plotDecorate = function( treeList, treeListClusters, query, cols=c( "white","red"), plotTree=TRUE, windowWidth=1000){

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

  # create new set of clusters including non-clusters
  clstComplete = rep(NA,nrow(fit[[1]]@correlation))
  names(clstComplete) = rownames(fit[[1]]@correlation)

  idx2 = match(names(clst), names(clstComplete))
  clstComplete[idx2] = clst

  if( any(!idx) & plotTree){
    cat("Cannot plot tree when features are dropped\n")
    plotTree = FALSE
  }

  # Make Views
  #############

  w = 0.6
  w2 = sqrt(2*w^2)
  # w2 = w

  # standard view
  heatmapVP <- viewport(width = unit(w2, "snpc"), height = unit(w2, "snpc"),
                        name="straight")

  # rotate correlation plot
  flipVP <- viewport(width = unit(w, "snpc"), height= unit(w, "snpc"), y=.6, angle=-45, name="flipVP")

  pushViewport(heatmapVP)

  # Title
  ########
  df = data.frame(query)
  main = paste0(df$seqnames, ':', format(df$start, big.mark=','), '-', format(df$end, big.mark=','))  
  title <- textGrob(main, 0.5, 1.05, gp=gpar(cex=1.0), name="title")

  # Plot Triangle
  ################

  cols = c("white","red")
  colsFun = colorRampPalette( (cols ))
  color = colsFun(100)

  flip=TRUE
  # Flip or not, determines way data is read into the display
  byrow<-ifelse(flip,FALSE,TRUE) #FALSE if flip=TRUE

  C = as.matrix(fit[[1]]@correlation )
  C[lower.tri(C, diag=TRUE)] = NA

  mybreak <- 0:length(color)/length(color)

  fill = as.character(cut(C^2,mybreak,labels=as.character(color), include.lowest=TRUE))

  ImageRect = makeImageRect( nrow(C), ncol(C), fill, name='heatmap', byrow=byrow, cutoff=windowWidth)   

  ImageRect <- editGrob(ImageRect, vp=flipVP)

  key = addLegend(rev(color), heatmapVP)
  
  heatMap <- gTree(children=gList(ImageRect, key), name="heatMap")

  # Plot segments
  ###############

  cols = rainbow( length(unique(clstComplete)))[clstComplete]

  N = length(clstComplete)

  figSegments = rectGrob(x=c(1:N)/N, y=.7, 
             width=1/N, height=1/20, 
             just=c("right", "top"), 
             gp=gpar(col=NA, fill=cols),
             name='d')

  # Combine and plot
  ###################

  if( plotTree ){

    Segs = gTree(children=gList(figSegments), name="segs")

    treeTmp = fit[[1]]@clust 

    # if tree is chac from adjclust, use correct() to get positive branch lengths
    if( is(treeTmp, 'chac') ){
      treeTmp = as.hclust(correct( treeTmp ))
    }

    figTree = ggdendrogram(treeTmp, labels=FALSE, leaf_labels=FALSE) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())

    treeVP <- viewport(y = unit(.8, "npc"), width = unit(.96, "snpc"), height = unit(w2/3, "snpc"), name="straight")

    plot.grob <- ggplotGrob(figTree)

    grid.newpage()
    grid.draw(gTree(children=gList(plot.grob, title), vp=treeVP))
    grid.draw(gTree(children=gList(figSegments), vp=heatmapVP))
    grid.draw(heatMap)
  }else{

    Segs = gTree(children=gList(figSegments, title), name="segs")

    grid.newpage()
    grid.draw(heatMap)
    grid.draw(gTree(children=gList(figSegments), vp=heatmapVP))
  }
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
#' # Choose cutoffs and return clusters
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







