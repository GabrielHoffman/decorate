
# library(ggplot2)
# n = 200
# p = 3
# data = matrix(rnorm(n*p), n, p)

# X = as.data.frame(data[1:100,])
# Y = as.data.frame(data[101:200,])

#' Scatter plot of all pairs of variables stratified by test variable
#'
#' Make a scatterplot for each pair of variables in X and Y.  Dataset is divided in two based on value in testVariable
#'
#' @param X data.frame of variables
#' @param Y data.frame of variables
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param size size of points
#' @param cols color to label samples in for two levels of testVariable
#' @param axisLabels either "show" to display axisLabels, "internal" for labels  in the diagonal plots, or "none" for no axis labels
#' @param title title
#' @param xlab xlab
#' @param ylab xlab
#'
#' @return ggplot2 of combined pairwise scatter plots
#'
#' @import ggplot2
#' @importFrom GGally ggmatrix
#' @importFrom utils combn
plotPairwiseScatter = function(X, Y, testVariable, size = 1, cols=c("darkgreen", "navy"), axisLabels = c("show", "internal", "none"), title=NULL, xlab = NULL, ylab = NULL){

    axisLabels <- match.arg( axisLabels )

    if( ! is.data.frame(X) ){
        stop("X must be of type data.frame")
    }    
    if( ! is.data.frame(Y) ){
        stop("Y must be of type data.frame")
    }
    if( ncol(X) != ncol(Y) ){
        stop("X and Y must have same number of columns")
    }
    if( ! is(testVariable, 'factor') ){
        stop("Entries in testVariable must be a factor with entries in two levels")
    }
    testVariable = droplevels(testVariable)
    if( nlevels(testVariable) != 2 ){
        stop("Entries in testVariable must be a factor with entries in two levels")
    }
    if( length(cols) != 2){
        stop("cols must have two entries")
    }

    # construct index pairs
    p = ncol(X)
    idx = t(combn(p, 2))
    idx = rbind(idx, idx[,2:1])
    idx = rbind(idx, cbind(seq_len(p), seq_len(p)))
    idx = idx[order(idx[,1], idx[,2]),]
    idx = data.frame(idx)
    idx$order = seq_len(nrow(idx))

    value = variable = a = b = NULL

    # plot all pairs
    ggpairsPlots = lapply( seq_len(nrow(idx)), function(k){

        if( idx[k,1] == idx[k,2]){
            df_x = data.frame(variable=levels(testVariable)[1], value = X[,idx[k,1]])
            df_y = data.frame(variable=levels(testVariable)[2], value = Y[,idx[k,1]])
            df = rbind(df_x, df_y)

            fig = ggplot(df, aes(value, color=variable)) + geom_density() + theme_bw(17) + scale_color_manual(values = cols) 
        }else if( idx[k,1] > idx[k,2]){

            df = data.frame(a = X[,idx[k,1]], b = X[,idx[k,2]])

            fig = ggplot(df, aes(a,b)) + geom_point(size=size, color=cols[1]) + theme_bw(17) + geom_abline(linetype=2) + geom_smooth(method='lm', se=FALSE, color="black")
            # + annotate("text", x=mean(df$a), y=mean(df$b), label=paste(k, paste(idx[k,], collapse=', '))) 
        }else{
            df = data.frame(a = Y[,idx[k,1]], b = Y[,idx[k,2]])

            fig = ggplot(df, aes(a,b)) + geom_point(size=size, color=cols[2]) + theme_bw(17) + geom_abline(linetype=2)  + geom_smooth(method='lm', se=FALSE, color="black")
        }
     fig
    })
    names(ggpairsPlots) = seq_len(length(ggpairsPlots))

    # reorder plots
    idx = idx[order(idx[,2], p-idx[,1]),]    

    columnLabels = colnames(X)
    labeller = "label_value"; switch = NULL; showStrips = TRUE;
    legend = NULL; cardinality_threshold = 15; progress = FALSE;

    fig = ggmatrix(plots = ggpairsPlots[idx$order], byrow = TRUE, nrow = p, ncol = p, xAxisLabels = (if (axisLabels ==
            "internal")
            NULL
        else rev(columnLabels)), yAxisLabels = (if (axisLabels ==
            "internal")
            NULL
        else columnLabels), labeller = labeller, switch = switch,
        showStrips = showStrips, showXAxisPlotLabels = identical(axisLabels,
            "show"), showYAxisPlotLabels = identical(axisLabels,
            "show"), title = title, xlab = xlab, ylab = ylab, legend=c(1,p))
    fig + theme(plot.title = element_text(hjust = 0.5)) 
}



#' Scatter plot of all pairs of variables stratified by test variable
#'
#' Scatter plot of all pairs of variables stratified by test variable
#'
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param peakIDs feature names to extract from rows of epiSignal
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param size size of points
#'
#' @return ggplot2 of combined pairwise scatter plots
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
#' treeListClusters = createClusters( treeList, method = "meanClusterSize", meanClusterSize=c(10, 30) )
#'
#' # get peak ID's from chr1, cluster 1
#' peakIDs = getFeaturesInCluster( treeListClusters, "chr20", 2, "30")
#'
#' # plot comparison of correlation matrices for peaks in peakIDs
#' #  where data is subset by metadata$Disease
#' plotScatterPairs( simData, peakIDs, metadata$Disease) + ggtitle("chr20: cluster 1")
#'
#' @import ggplot2
#' @importFrom GGally ggmatrix
#' @export
plotScatterPairs = function(epiSignal, peakIDs, testVariable, size=1){

  if( length(testVariable) != ncol(epiSignal) ){
    stop("Number of columns in epiSignal must equal number of entries in testVariable")
  }
  if( ! is(testVariable, 'factor') ){
    stop("Entries in testVariable must be a factor with entries in two levels")
  }
  testVariable = droplevels(testVariable)
  if( nlevels(testVariable) != 2 ){
    stop("Entries in testVariable must be a factor with entries in two levels")
  }

  # Divide data into two sets
  set1 = which(testVariable == levels(testVariable)[1])
  set2 = which(testVariable == levels(testVariable)[2])

  # get features
  Var1 = Var2 = value = NA
  Y1 = t(epiSignal[peakIDs, set1])
  Y2 = t(epiSignal[peakIDs, set2])

  Y1 = data.frame( Y1 )
  Y2 = data.frame( Y2 )

  plotPairwiseScatter( Y1, Y2, testVariable, size=size, axisLabels='show')
}


# source("/Users/gabrielhoffman/workspace/repos/decorate/R/plotCorrPairs.R");  plotScatterPairs( simData, peakIDs[1:4], metadata$Disease) 

# dev.new()
# plot(t(simData[peakIDs[1:2],  metadata$Disease==1]), pch=20)

# plotCompareCorr( simData, peakIDs[1:4], metadata$Disease) + ggtitle("chr1: cluster 1")




