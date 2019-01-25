


#' Evaluate Differential Correlation
#' 
#' Evaluate Differential Correlation between two subsets of data
#' 
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param gr GenomciRanges corresponding to the rows of epiSignal
#' @param clustList list of cluster assignments
#' @param npermute array of two entries with min and max number of permutations
# @param BPPARAM parameters for parallel evaluation by chromosome
#' @param mc.cores number of cores to permutations
#' 
#' @return list of result by chromosome and clustList
#'
#' @details
#' Correlation sturucture between two subsets of the data is evaluated with sparse-Leading-Eigenvalue-Driven (sLED) test:
#' Zhu, Lingxue, Jing Lei, Bernie Devlin, and Kathryn Roeder. 2017. Testing high-dimensional covariance matrices, with application to detecting schizophrenia risk genes. Annals of Applied Statistics. 11:3 1810-1831. doi:10.1214/17-AOAS1062
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
#' # Simulate variable to split dataset by
#' set.seed(1)
#' metadata = data.frame( Disease = factor(sample(0:1, ncol(simData), replace=TRUE)))
#'
#' # Evaluate Differential Correlation between two subsets of data
#' sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute=20)
#'
#' # get summary of results
#' df = summary( sledRes )
#'
#' # print results
#' head(df)
#' 
#' @importFrom GenomicRanges GRanges
#' @import BiocParallel
#' @export
#' @docType methods
#' @rdname evalDiffCorr-methods
setGeneric("evalDiffCorr", function(epiSignal, testVariable, gr, clustList , npermute = c(100, 10000), mc.cores=4) standardGeneric("evalDiffCorr"))
	# BPPARAM = SerialParam()) standardGeneric("evalDiffCorr"))

#' @import limma
#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,EList,ANY,GRanges,list,numeric,ANY-method
setMethod("evalDiffCorr", c("EList", "ANY", "GRanges", "list", "ANY", "ANY"), 
	function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000), mc.cores=4){#, BPPARAM = SerialParam()){
		.evalDiffCorr( epiSignal$E, testVariable, gr, clustList, npermute, mc.cores) #BPPARAM )
	})


#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,matrix,ANY,GRanges,list,numeric,ANY-method
setMethod("evalDiffCorr", c("matrix", "ANY", "GRanges", "list", "ANY", "ANY"), 
	function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000), mc.cores=4){#, BPPARAM = SerialParam()){
		.evalDiffCorr( epiSignal, testVariable, gr, clustList, npermute, mc.cores) #BPPARAM )
	})


#' An S4 class that stores results of sLED analysis
#' @slot list
#' @export
setClass("sLEDresults", representation("list"))

#' @import sLED
.evalDiffCorr = function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000), mc.cores=4){#, BPPARAM = SerialParam()){

	if( nrow(epiSignal) != length(gr)){
		stop("Number of rows in epiSignal must equal number of entries in gr")
	}
	if( length(testVariable) != ncol(epiSignal) ){
		stop("Number of columns in epiSignal must equal number of entries in testVariable")
	}
	if( ! is(testVariable, 'factor') || nlevels(testVariable) != 2 ){
		stop("Entries in testVariable must be a factor with two levels")
	}
	if( min(table(testVariable)) < 20){
		stop("Need at last 20 samples in smallest class of testVariable")
	}
	if( length(npermute) != 2 || npermute[1] >= npermute[2] ){
		stop("npermute must have two entries: min and max permutations")
	}

	allClusters = unlist(lapply(names(clustList), function(x) paste0(x, '_', unique(clustList[[x]]))))

	set1 = which(testVariable == levels(testVariable)[1])
	set2 = which(testVariable == levels(testVariable)[2])

	# for each chromosome
	res = lapply( names(clustList), function(chrom){

		chromClust = clustList[[chrom]]
	  	clustValTab = table(chromClust)
	  	clustValTab = clustValTab[clustValTab >3]
	  	clustIDLst = names(clustValTab)

	  # for each cluster in the chromosome
	  sledRes = lapply( clustIDLst, function( clstLabel ){

	    i = match(paste0(chrom, '_', clstLabel), allClusters)

		cat('\r', i, ' / ', length(allClusters))

	    # get names of peaks in this cluster: clstLabel
	    peakIDs = names(chromClust)[chromClust == clstLabel]

	    if( length(peakIDs) > 3 ){

	    	# get two subsets of data
	      	Y1 = t(epiSignal[peakIDs,set1])
	      	Y2 = t(epiSignal[peakIDs,set2])

	      	# perform permutations until p-values is precise enough
	      	# if not precise enough
	      	a = log10(npermute[1])
	      	b = log10(npermute[2])
	      	permArray = 10^seq(a,b, length.out=b-a +1)

	      	for( nperm in permArray){
		      	# compare correlation structure with sLED
		      	res = sLED(X=Y1, Y=Y2, npermute=nperm, verbose=FALSE, mc.cores=mc.cores, useMC=mc.cores>1)

		      	if( res$pVal * nperm > 10){
		      		break
		      	}
		    }

	    }else{
	      	res = NULL
	    }

	    if( i %% 100 == 0){
			gc()
		}

	    res
	  })
	  names(sledRes) = clustIDLst
	  sledRes
	})
	names(res) = names(clustList)
	
	cat("\n")
	new("sLEDresults", res)
}

#' Summarize sLED analysis
#' 
#' extract statistic and p-value for each cluster
#' 
#' @param object sLEDresults
#' 
#' @return data.frame
#'
#' @export
setMethod("summary", "sLEDresults", function( object ){

	res = lapply( names(object), function(chrom){
		res = lapply( names(object[[chrom]]), function(clst){
			df = object[[chrom]][[clst]]
			data.frame(chrom=chrom, cluster=clst, pValue = df$pVal, stat=df$stat, n.pem=length(df$Tn.permute))
		})
		do.call("rbind", res)
	})
	res = do.call("rbind", res)

	res$p.adjust = p.adjust( res$pValue, "fdr" )
	res
})





