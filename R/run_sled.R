


#' Evaluate Differential Correlation
#' 
#' Evaluate Differential Correlation between two subsets of data
#' 
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param gr GenomciRanges corresponding to the rows of epiSignal
#' @param clustList list of cluster assignments
#' @param npermute array of two entries with min and max number of permutations
#' @param adj.beta parameter for sLED
#' @param BPPARAM parameters for parallel evaluation by chromosome
#' 
#' @return list of result by chromosome and clustList
#'
#' @details
#' Correlation sturucture between two subsets of the data is evaluated with sparse-Leading-Eigenvalue-Driven (sLED) test:
#'
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
#' sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute=c(20, 200))
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
setGeneric("evalDiffCorr", function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000), adj.beta=-1, BPPARAM = SerialParam()) standardGeneric("evalDiffCorr"))

#' @import limma
#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,EList,ANY,GRanges,list,ANY,ANY,ANY-method
setMethod("evalDiffCorr", c("EList", "ANY", "GRanges", "list", "ANY", 'ANY', "ANY"), 
	function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000),  adj.beta=-1, BPPARAM = SerialParam()){
		.evalDiffCorr( epiSignal$E, testVariable, gr, clustList, npermute, adj.beta, BPPARAM)
	})


#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,matrix,ANY,GRanges,list,ANY,ANY,ANY-method
setMethod("evalDiffCorr", c("matrix", "ANY", "GRanges", "list", "ANY", 'ANY',"ANY"), 
	function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000), adj.beta=-1, BPPARAM = SerialParam()){
		.evalDiffCorr( epiSignal, testVariable, gr, clustList, npermute, adj.beta, BPPARAM)
	})


#' An S4 class that stores results of sLED analysis
#' @slot list
#' @export
setClass("sLEDresults", representation("list"))

#' @import sLED
runSled = function(i, dfClustUnique, dfClust, epiSignal, set1, set2, npermute){

	CHROM = dfClustUnique$chrom[i]
	CLST = as.character(dfClustUnique$cluster[i])

	cat(CHROM, CLST, '\n')

	peakIDs = dfClust[(chrom==CHROM) & (cluster==CLST),peak]

	if( length(peakIDs) > 3 ){

		# get two subsets of data
	  	Y1 = t(epiSignal[peakIDs,set1])
	  	Y2 = t(epiSignal[peakIDs,set2])

	  	# perform permutations until p-values is precise enough
	  	# if not precise enough
	  	a = log10(npermute[1])
	  	b = log10(npermute[2])
	  	permArray = 10^seq(a,b, length.out=b-a +1)

	  	for( nperm in round(permArray) ){
	      	# compare correlation structure with sLED
	      	Y1_scale = scale(Y1)
	      	Y2_scale = scale(Y2)

	      	res = sLED(X=Y1_scale, Y=Y2_scale, npermute=nperm, verbose=FALSE, mc.cores=1, useMC=FALSE, adj.beta=adj.beta)

	      	if( res$pVal * nperm > 10){
	      		break
	      	}
	    }

	}else{
	  	res = list(pVal=NA, stats=NA)
	}
	res$chrom = CHROM
	res$cluster = CLST
	res
}



icount2 = function (count){
    if (missing(count))
        count <- NULL
    else if (!is.numeric(count) || length(count) != 1)
        stop("count must be a numeric value")
    i <- 0L
    nextEl <- function(){
    	if( is.null(i) )    		
        	(i <<- NULL)
        else if (is.null(count) || i < count)
            (i <<- i + 1L)
        else 
        	(i <<- NULL)
    }
    it <- list(nextElem = nextEl)
    class(it) <- c("abstractiter", "iter")
    it
}

#' @import iterators
#' @importFrom data.table data.table
clustIter = function( dfClustUnique, dfClust, epiSignal, set1, set2 ){

	n_clusters = nrow( dfClustUnique )

	xit = icount2( n_clusters )

	nextEl <- function( ) {

		i <- nextElem( xit )

		if( is.null(i) || i > n_clusters){
			res = NULL
		}else{
			CHROM = dfClustUnique$chrom[i]
			CLST = dfClustUnique$cluster[i]

			peakIDs = dfClust[(chrom==CHROM) & (cluster==CLST),peak]

			# get two subsets of data
		  	Y1 = t(epiSignal[peakIDs,set1,drop=FALSE])
		  	Y2 = t(epiSignal[peakIDs,set2,drop=FALSE])

		  	if( ncol(Y1) != ncol(Y2) ){
		  		stop("ncol(Y1) != ncol(Y2): ", ncol(Y1), ncol(Y2), "\n", i)
		  	}

		  	res = list(Y1=Y1, Y2=Y2, CHROM=CHROM, CLST=CLST, i=i)
		}
		res
	}

	it <- list(nextElem = nextEl)
    class(it) <- c(  "abstractiter", "iter") 
	it
}

#' @import sLED
runSled2 = function( itObj, npermute, adj.beta){

	ncol1 = ncol(itObj$Y1)
	ncol2 = ncol(itObj$Y2)

	if( min(ncol1, ncol2) > 3 ){

	  	# perform permutations until p-values is precise enough
	  	# if not precise enough
	  	a = log10(npermute[1])
	  	b = log10(npermute[2])
	  	permArray = 10^seq(a,b, length.out=b-a +1)

	  	for( nperm in round(permArray) ){
	      	# compare correlation structure with sLED
	      	res = sLED(X=scale(itObj$Y1), Y=scale(itObj$Y2), npermute=nperm, verbose=FALSE, mc.cores=1, useMC=FALSE, adj.beta=adj.beta)

	      	if( res$pVal * nperm > 10){
	      		break
	      	}
	    }

	}else{
	  	res = list(pVal=NA, stats=NA)
	}
	res$chrom = itObj$CHROM
	res$cluster = itObj$CLST
	res
}



#' @import BiocParallel
#' @importFrom data.table data.table
.evalDiffCorr = function(epiSignal, testVariable, gr, clustList, npermute = c(100, 10000), adj.beta=-1, BPPARAM = SerialParam()){

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
	if( any(!names(treeListClusters) %in% seqnames(peakLocs2)@values) ){
		stop("Chromosomes from treeListClusters must be in peakLocs2")
	} 


	allClusters = unlist(lapply(names(clustList), function(x) paste0(x, '_', unique(clustList[[x]]))))

	# Divide data into two sets
	set1 = which(testVariable == levels(testVariable)[1])
	set2 = which(testVariable == levels(testVariable)[2])

	# create single data.table with cluster and peak info
	dfClust = lapply( names(clustList), function(chrom){
		data.frame(chrom=chrom, cluster=clustList[[chrom]], peak=names(clustList[[chrom]]), stringsAsFactors=FALSE)
	})
	dfClust = data.table(do.call('rbind', dfClust))

	# only get peaks that are in the the epiSignal dataset
	dfClust = dfClust[peak %in% rownames(epiSignal),]

	# get unique chrom and 
	dfClustUnique = unique(dfClust[,1:2])

	cat("# Clusters:", nrow(dfClustUnique), '\n')

	# Evaluate statistics with permutations
	#######################################

	# res = mclapply( seq_len(nrow(dfClustUnique)), runSled, dfClustUnique, dfClust, epiSignal, set1, set2, npermute=c(10, 100), mc.cores=1)

	# lapply(1, runSled, dfClustUnique, dfClust, epiSignal, set1, set2, npermute)

	 # runSled(1, dfClustUnique, dfClust, epiSignal, set1, set2, npermute)

	# combinedResults = bplapply( seq_len(nrow(dfClustUnique)), runSled, dfClustUnique, dfClust, epiSignal, set1, set2, npermute, BPPARAM=BPPARAM)
	
	# run with iterators
	it = clustIter( dfClustUnique, dfClust, epiSignal, set1, set2  )
	
	combinedResults = bpiterate( it$nextElem, runSled2, npermute, adj.beta, BPPARAM=BPPARAM)

	# return list of lists
	######################
	resList = list()
	chromArray = unique(vapply(combinedResults, function(x) x$chrom, "character"))
	for( chrom in chromArray){

		id = vapply(combinedResults, function(x) x$cluster, numeric(1))

		idx = vapply(combinedResults, function(x){ x$cluster %in% id && x$chrom == chrom}, logical(1)) 

		resList[[chrom]] = combinedResults[idx]
		names(resList[[chrom]]) = vapply(resList[[chrom]], function(x) x$cluster, numeric(1))

	}
	names(resList) = chromArray
	
	cat("\n")
	new("sLEDresults", resList)
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
			data.frame(chrom=chrom, cluster=clst, pValue = df$pVal, stat=df$stat, n.perm=length(df$Tn.permute))
		})
		do.call("rbind", res)
	})
	res = do.call("rbind", res)

	res$p.adjust = p.adjust( res$pValue, "fdr" )
	res[order(res$pValue),]
})





