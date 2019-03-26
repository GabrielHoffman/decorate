
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


#' Convert correlation matrix into triangle plot
#'
#' Adapted from LDheatmap
#'
#' @param nrow nrow(C)
#' @param ncol ncol(C)
#' @param cols entries of C converted to color
#' @param name name of the plot
#' @param byrow process C by row
# @param cutoff retain correlations up to cutoff features away 
#'
#' @return rectGrob 
makeImageRect <- function(nrow, ncol, cols, name, byrow=TRUE) {
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

  # fast plot
  # i = (right-top)^2 < cutoff/nrow

  # rectGrob(x=right[i], y=top[i], 
  #          width=1/ncol, height=1/nrow, 
  #          just=c("right", "top"), 
  #          gp=gpar(col=NA, fill=cols[i], alpha=1),
  #          name=name)

  rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols, alpha=1),
           name=name)
}


addLegend <- function(color, vp){
  ImageRect<- makeImageRect(2,length(color), cols=c(rep(NA,length(color)),color[length(color):1]), "colorKey")

  keyVP <- viewport(x=.94, y=.03, height=.05, width=.2, just=c("right","bottom"), name="keyVP")
  
  title<-textGrob("abs(Correlation)", x=0.5, y=1.25, name="title", gp=gpar(cex=0.8))
  
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
#' @param ensdb ENSEMBL database object like EnsDb.Hsapiens.v86 
#' @param treeList hierarchical clustering of each chromosome from runOrderedClusteringGenome()
#' @param treeListClusters assign regions to clusters after cutting tree with createClusters()
#' @param featurePositions GRanges object storing location of each feature
#' @param query GRanges object indiecating region to plot
#' @param cols vector of two colors
#' @param showTree show tree from hierachical clustering
#' @param showGenes plot genes
#' @param splice_variants if TRUE, show multiple transcripts from the same gene
#' @param non_coding if TRUE, also show non-coding genes
#'
#' @return ggplot2 of cluster assignments and correlation between peaks
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
#' treeListClusters = createClusters( treeList )
#' 
#' # Plot correlations and clusters in region defind by query
#' query = range(simLocation)
#' 
#' plotDecorate( ensdb, treeList, treeListClusters, simLocation, query)
#'
#' @import ggplot2
#' @import grid
#' @import Biobase
#' @import GenomicRanges
#' @importFrom labeling extended
#' @importFrom stats as.hclust
#' @importFrom ggdendro ggdendrogram
#' @importFrom reshape2 melt
#' @importFrom cowplot ggdraw draw_grob
#' @importFrom adjclust correct
#' @importFrom grDevices rainbow colorRampPalette
#' @export
plotDecorate = function( ensdb, treeList, treeListClusters, featurePositions, query, cols=c( "lightyellow","red"), showTree=FALSE, showGenes=TRUE, splice_variants=FALSE, non_coding=FALSE){

  if( ! is(query, "GRanges")){
    stop("query must be GRanges object")
  }

  if( ! is(featurePositions, "GRanges")){
    stop("featurePositions must be a GRanges object")
  }

  if( is.null(names(featurePositions)) ){
    stop("featurePositions must have identifiers for each interval accessable using names(featurePositions)")
  }

  if( ! is(ensdb, "EnsDb")){
    stop("ensdb must be EnsDb object")
  }

  if( length(query) > 1){
    stop("Can only query one interval")
  }

  # drop empty chromosomes
  chrNameCount = table(seqnames(query))

  query = dropSeqlevels( query, names(chrNameCount[chrNameCount==0]))

  # get clusters in query
  fit = getSubset( treeList, query )

  if( length(fit) == 0){
    stop("No data left after filtering by query region")
  }

  if( length(fit[[1]]@correlation) < 2 ){
    stop("Must have at least two features in this region")
  }

  if( ! is(treeListClusters, 'list') && ! is(treeListClusters[[1]], 'epiclustDiscreteList') ){
    stop("treeListClusters must be the result of createClusters") 
  }

  # check that features labels are also found in featurePositions
  idx = match( fit[[1]]@clust$labels, names(featurePositions))

  if( any(is.na(idx)) ){
    stop("There are ", format(sum(is.na(idx)), big.mark=','), " features in treeListClusters that are not found in featurePositions")
  }

  chrom = names(fit)

  # if( showTree){
  #   cat("Dropping tree plot tree unless entire chromosome is plotted\n")
  #   showTree = FALSE
  # }

  ##############
  # Plot Genes #
  ##############

  if( showGenes ){

    if( width(query) > 5e7){
      stop("Range is too wide (>5e7 bp) to plot genes: ", width(query), " bp.\nYou can plot this range with showGenes=FALSE")
    }

    # library(LDheatmap)
    # transcripts <- plotGenes(start(query), end(query), seqnames(query),
        # genome="hg38", plot_lines_distance = 0.04, splice_variants=FALSE, non_coding=FALSE)

    transcripts <- plotEnsGenes(ensdb, start(query), end(query), seqnames(query), plot_lines_distance = 0.03, non_coding=non_coding, splice_variants=splice_variants)

    startHeight = attr(transcripts, 'height') 
  }else{
    startHeight = 0.01
  }

  #################
  # Plot segments #
  #################

  y_feature_locs = 1 - startHeight
  feat_mark_y = 0.07
  cluster_mark_y = 0.03

  y_pos = y_feature_locs - 0.02

  n_features = length(fit[[1]]@clust$labels)

  clustColsLst = lapply(treeListClusters, function(tlc){
     clst = tlc[[chrom]]
    idx = names(clst) %in% fit[[1]]@clust$labels
    clst = clst[idx]

    # create new set of clusters including non-clusters
    clstComplete = rep(NA,nrow(fit[[1]]@correlation))
    names(clstComplete) = rownames(fit[[1]]@correlation)

    idx2 = match(names(clst), names(clstComplete))
    clstComplete[idx2] = clst
    rainbow( length(unique(clstComplete)))[as.factor(clstComplete)]
  })

  # remove lines with no segments
  clustColsLst = lapply(clustColsLst, function(x){
    ret = NA
    if(any(!is.na(x))){
      ret = x
    }
    ret
  })
  clustColsLst = clustColsLst[lengths(clustColsLst)>1]

  N = length(clustColsLst)
  xval = seq(0, 1, length.out=n_features+1)[1:n_features]
  incr = 0.013 
  # yval = sort(rep(rep(y_feature_locs-cluster_mark_y-feat_mark_y+incr*0:(N-1)), n_features))
  yval = sort(rep(rep(y_pos - 1:N*incr), n_features))
  segCols = unlist(clustColsLst)

  if( length(xval) > 0){
    figSegments = segmentsGrob( x0=rep(xval,N), 
          x1 = rep(xval+1/n_features,N), y0=yval, y1=yval, gp=gpar(lwd=3, col=segCols, lineend="butt"))
  }else{
    figSegments = segmentsGrob()
  }

  ####################################
  # Plot feature locations  and grid # 
  ####################################

  # subset of featurePositions based on query
  fnd = findOverlaps(featurePositions, query)
  featurePositions = featurePositions[fnd@from]

  pos1 = (start(featurePositions) - start(query)) / width(query) 
  pos2 = (end(featurePositions) - start(query)) / width(query)
  midpoint = pos1 + (pos2 -pos1)/2

  # cbind(pos1, pos2, pos1 + (pos2-pos1)/2, midpoint)

  figFeatLocations = segmentsGrob( x0=pos1, x1=pos2, y0=y_pos, y1=y_pos, gp=gpar(lwd=3, col='navy', lineend="butt"))

  xval = seq(0, 1, length.out=n_features+1)[1:n_features] + 1/n_features/2

  figPos = segmentsGrob( x0=xval, x1=midpoint, y0=max(yval), y1=y_pos, gp=gpar(lwd=1, col="grey80", lineend="butt"))

  # Make Views
  #############

  ylocation = min(yval) * .92

  w = 0.6
  w2 = sqrt(2*w^2)
  # w2 = w

  # standard view
  heatmapVP <- viewport(width = unit(w2, "snpc"), height = unit(w2, "snpc"), name="straight")

  # rotate correlation plot
  flipVP <- viewport(width = unit(w, "snpc"), height= unit(w, "snpc"), y=ylocation, angle=-45, name="flipVP", gp=gpar(fill="black"))

   flipVP2 <- viewport(width = unit(0.705, "snpc"), height= unit(0.705, "snpc"), angle=-45, name="flipVP2", y=ylocation)

  # Title
  ########
  df = data.frame(query)
  main = paste0(df$seqnames, ':', format(df$start, big.mark=','), '-', format(df$end, big.mark=','))  
  title <- textGrob(main, 0.5, 1.05, gp=gpar(cex=1.0), name="title")

  # Plot Triangle
  ################

  colsFun = colorRampPalette( (cols ))
  color = c('white', colsFun(1000))

  flip=TRUE
  # Flip or not, determines way data is read into the display
  byrow<-ifelse(flip,FALSE,TRUE) #FALSE if flip=TRUE

  C = as.matrix(fit[[1]]@correlation )
  C[lower.tri(C, diag=TRUE)] = NA

  n_features = nrow(C)

  if( !identical(rownames(C), names(featurePositions)) ){
    cat("Processing subset of features in query...\n")
    showTree = FALSE
  }

  mybreak <- 0:length(color)/length(color)
  mybreak[2] = 1e-7

  fill = as.character(cut(abs(C),mybreak,labels=as.character(color), include.lowest=TRUE))

  ImageRect = makeImageRect( nrow(C), ncol(C), fill, name='heatmap', byrow=byrow)
  ImageRect <- editGrob(ImageRect, vp=flipVP)

  key = addLegend(rev(color), heatmapVP)
  
  heatMap <- gTree(children=gList(ImageRect, key), name="heatMap")
  
  ##############
  # Plot title #
  ##############

  txt = with(data.frame(query), paste0(seqnames, ':',format(start,big.mark=','),'-',format(end, big.mark=',')))

  plot_title <- textGrob(txt,
        gp = gpar(fontsize = 12, fontfamily = "mono"), y = 1,
        just = c("centre", "center"), name = "gene_plot_title",
        default.units = "native")

  # Combine and plot
  ###################

  if( showTree ){
    Segs = gTree(children=gList(figSegments), name="segs")

    treeTmp = fit[[1]]@clust 

    # if tree is chac from adjclust, use correct() to get positive branch lengths
    if( is(treeTmp, 'chac') ){
      treeTmp = as.hclust(correct( treeTmp ))
    }

    figTree = ggdendrogram(treeTmp, labels=FALSE, leaf_labels=FALSE, size=.2) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), plot.background = element_rect(fill="white", color = NA))

    treeVP <- viewport(y = unit(.8, "npc"), width = unit(.96, "snpc"), height = unit(w2/3, "snpc"), name="straight")

    plot.grob <- ggplotGrob(figTree)

    grid.newpage()
    ggdraw(heatMap) + theme(plot.background = element_rect(fill="white", color = NA)) + 
      draw_grob(gTree(children=gList(plot.grob, title), vp=treeVP)) + 
      draw_grob(gTree(children=gList(figSegments), vp=heatmapVP))  
      # draw_grob(gTree(children=gList(plot_title), vp=heatmapVP))
  }else{

    # grid.newpage()
    # fig = ggdraw()

    if( showGenes ){
      # fig = fig + draw_grob(gTree(children=gList(transcripts), vp=heatmapVP))
      fig = ggdraw(heatMap) + draw_grob(gTree(children=gList(transcripts, figPos, figFeatLocations, figSegments, plot_title), vp=heatmapVP)) 
    }else{

     fig =  ggdraw(heatMap) + draw_grob(gTree(children=gList(figPos, figFeatLocations, figSegments, plot_title), vp=heatmapVP)) 
    }

    # fig + draw_grob(gTree(children=gList(figPos), vp=heatmapVP)) +
    #   draw_grob(gTree(children=gList(figFeatLocations), vp=heatmapVP)) + 
    #   draw_grob(gTree(children=gList(figSegments), vp=heatmapVP))  + 
    #   draw_grob(heatMap) + draw_grob(gTree(children=gList(plot_title), vp=heatmapVP))
    fig
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
#' treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c( 10, 20) )
#'
#' # Simulate variable to split dataset by
#' set.seed(1)
#' metadata = data.frame( Disease = factor(sample(0:1, ncol(simData), replace=TRUE)))
#'
#' # get peak ID's from chr1, cluster 1
#' peakIDs = getFeaturesInCluster( treeListClusters, "chr1", 1, "10")
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

  ggplot(df, aes(Var1, Var2)) + geom_tile(aes(color=value, fill=value)) + scale_color_gradientn(name = "Correlation", colours=cols, limits=c(-1,1), na.value="grey") + scale_fill_gradientn(name = "Correlation", colours=cols, limits=c(-1,1), na.value="grey")  + theme_void() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="bottom") + annotate(geom="text", x=c(1-.2,N+.2), y=c(N,1), label=levels(testVariable), size=size, hjust=c(0,1), vjust=c(.5,.5))
}


