


#' Evaluate Differential Correlation
#' 
#' Evaluate Differential Correlation between two subsets of data
#' 
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param gRanges GenomciRanges corresponding to the rows of epiSignal
#' @param clustList list of cluster assignments
#' @param npermute array of two entries with min and max number of permutations
#' @param adj.beta parameter for sLED
#' @param rho a large positive constant such that A(X)-A(Y)+diag(rep(rho,p)) is positive definite. Where p is the number of features
#' @param sumabs.seq sparsity parameter
#' @param BPPARAM parameters for parallel evaluation
#' @param method "sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob"
#' @param method2 Specify type of correlation: "pearson", "kendall", "spearman"
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
#' # Plot correlations and clusters in region defind by query
#' query = range(simLocation)
#' 
#' # Plot clusters
#' plotDecorate( ensdb, treeList, treeListClusters, simLocation, query)
#'
#' # Evaluate Differential Correlation between two subsets of data
#' sledRes = evalDiffCorr( simData, metadata$Disease, simLocation, treeListClusters, npermute=c(20, 200, 2000))
#'
#' # get summary of results
#' df = summary( sledRes )
#'
#' # print results
#' head(df)
#'
#' # extract peak ID's from most significant cluster
#' peakIDs = getFeaturesInCluster( treeListClusters, df$chrom[1], df$cluster[1], "20")
#'
#' # plot comparison of correlation matrices for peaks in peakIDs
#' #  where data is subset by metadata$Disease
#' main = paste0(df$chrom[1], ': cluster ', df$cluster[1])
#' plotCompareCorr( simData, peakIDs, metadata$Disease) + ggtitle(main)
#' 
#' @importFrom GenomicRanges GRanges
#' @import BiocParallel
#' @export
#' @docType methods
#' @rdname evalDiffCorr-methods
setGeneric("evalDiffCorr", function(epiSignal, testVariable, gRanges, clustList, npermute = c(100, 10000, 100000), adj.beta=0, rho = 0, sumabs.seq = 1, BPPARAM = SerialParam(), method=c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob" ), method2=c("pearson", "kendall", "spearman")) standardGeneric("evalDiffCorr"))

#' @import limma
#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,EList,ANY,GRanges,list,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method
setMethod("evalDiffCorr", c("EList", "ANY", "GRanges", "list", "ANY", "ANY", 'ANY', "ANY", "ANY", "ANY", "ANY"), 
	function(epiSignal, testVariable, gRanges, clustList, npermute = c(100, 10000, 100000), adj.beta=0, rho = 0, sumabs.seq = 1, BPPARAM = SerialParam(), method=c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob" ), method2=c("pearson", "kendall", "spearman")){
		.evalDiffCorr( epiSignal$E, testVariable, gRanges, clustList, npermute, adj.beta, rho, sumabs.seq, BPPARAM, method, method2)
	})


#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,matrix,ANY,GRanges,list,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method
setMethod("evalDiffCorr", c("matrix", "ANY", "GRanges", "list", "ANY", "ANY", 'ANY', "ANY", "ANY", "ANY", "ANY"), 
	function(epiSignal, testVariable, gRanges, clustList, npermute = c(100, 10000, 100000), adj.beta=0, rho = 0, sumabs.seq = 1, BPPARAM = SerialParam(), method=c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob" ), method2=c("pearson", "kendall", "spearman")){
		.evalDiffCorr( epiSignal, testVariable, gRanges, clustList, npermute, adj.beta, rho, sumabs.seq, BPPARAM, method, method2)
	})


#' @import BiocParallel
#' @export
#' @rdname evalDiffCorr-methods
#' @aliases evalDiffCorr,data.frame,ANY,GRanges,list,ANY,ANY,ANY,ANY,ANY,ANY,ANY-method
setMethod("evalDiffCorr", c("data.frame", "ANY", "GRanges", "list", "ANY", "ANY", 'ANY', "ANY", "ANY", "ANY", "ANY"), 
	function(epiSignal, testVariable, gRanges, clustList, npermute = c(100, 10000, 100000), adj.beta=0, rho = 0, sumabs.seq = 1, BPPARAM = SerialParam(), method=c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob" ), method2=c("pearson", "kendall", "spearman")){
		.evalDiffCorr( epiSignal, testVariable, gRanges, clustList, npermute, adj.beta, rho, sumabs.seq, BPPARAM, method, method2)
	})

#' An S4 class that stores results of sLED analysis
#'
#' @slot .Data list of sLED results
#' @export
setClass("sLEDresults", representation("list"))

#' Summarize sLED analysis
#' 
#' extract statistic and p-value for each cluster
#' 
#' @param object sLEDresults
#' 
#' @return data.frame
#'
#' @importFrom stats p.adjust
#' @export
setMethod("summary", "sLEDresults", function( object ){

	res = lapply( names(object), function(chrom){
		res = lapply( names(object[[chrom]]), function(clstKey){
			df = object[[chrom]][[clstKey]]
			if( is.null(df$sign) || is.na(df$sign) ){
				df$sign = "0"
			}
			id = unlist(strsplit(clstKey,' '))[1]
			clust = unlist(strsplit(clstKey,' '))[2]
			# multiply by direction of effect

			sgn = switch( df$sign, "pos" = 1, "0" = 0, "neg" = -1)
			data.frame(id=id, chrom=chrom, cluster=clust, pValue = df$pVal, stat=df$stats * sgn, n.perm=length(df$Tn.permute), stringsAsFactors=FALSE)
		})
		do.call("rbind", res)
	})
	res = do.call("rbind", res)

	res$p.adjust = p.adjust( res$pValue, "fdr" )
	res = res[order(res$pValue),]
	rownames(res) = c()
	res
})


setMethod("print", "sLEDresults", function( x ){
	cat("Hypothesis tests of each cluster\n\n")
	for( chrom in names(x) ){
	    nTest = length(x[[chrom]])
	    cat(paste0(chrom, ': Tests of'),nTest, "clusters\n" )
  	}
  	cat('\n')
})



setMethod("show", "sLEDresults", function( object ){
	print( object )
})



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
clustIter = function( dfClustUnique, dfClust, epiSignal, testVariable ){

	n_clusters = nrow( dfClustUnique )

	xit = icount2( n_clusters )

	nextEl <- function( ) {

		i <- nextElem( xit )

		if( is.null(i) || i > n_clusters){
			res = NULL
		}else{
			ID = dfClustUnique$id[i]
			CHROM = dfClustUnique$chrom[i]
			CLST = dfClustUnique$cluster[i]

			id = chrom = cluster = peak = NA
			peakIDs = dfClust[(id==ID) & (chrom==CHROM) & (cluster==CLST),peak]

			# get two subsets of data
		  	# Y1 = t(epiSignal[peakIDs,set1,drop=FALSE])
		  	# Y2 = t(epiSignal[peakIDs,set2,drop=FALSE])

		  	# if( ncol(Y1) != ncol(Y2) ){
		  	# 	stop("ncol(Y1) != ncol(Y2): ", ncol(Y1), ncol(Y2), "\n", i)
		  	# }

		  	res = list(Y 	= t(epiSignal[peakIDs,,drop=FALSE]),
		  			group 	= testVariable, 
		  			ID 		= ID, 
		  			CHROM 	= CHROM, 
		  			CLST 	= CLST, 
		  			i 		= i)
		}
		res
	}

	it <- list(nextElem = nextEl)
    class(it) <- c(  "abstractiter", "iter") 
	it
}

#' @import sLED
runSled2 = function( itObj, npermute, adj.beta, rho, sumabs.seq, BPPARAM){

	Y1 = itObj$Y[itObj$group==1,]
	Y2 = itObj$Y[itObj$group==1,]

	if( ncol(itObj$Y) >= 3 ){

	  	# perform permutations until p-values is precise enough
	  	# if not precise enough
	  	a = log10(npermute[1])
	  	b = log10(npermute[2])
	  	permArray = 10^seq(a,b, length.out=b-a +1)

	  	for( nperm in round(permArray) ){
	      	# compare correlation structure with sLED
	      	res = .sLED(X=Y1, Y=Y2, npermute=nperm, verbose=FALSE, adj.beta=adj.beta, rho=rho, sumabs.seq=sumabs.seq, BPPARAM=BPPARAM)

	      	if( res$pVal * nperm > 10){
	      		break
	      	}
	    }

	}else{
	  	res = list(pVal=NA, stats=NA, count=NA)
	}
	res$id = itObj$ID
	res$chrom = itObj$CHROM
	res$cluster = itObj$CLST
	res
}




#' Test difference between two correlation matricies
#'
#' Test difference between two correlation matricies using one of 5 tests
#'
# @param C1 correlation matrix
# @param C2 correlation matrix
# @param N1 number of samples use to estimate C1
# @param N2 number of samples use to estimate C2
#' @param Y data matrix
#' @param group a factor defining groups
#' @param method Specify test: "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney" "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob"
#' @param method2 Specify type of correlation: "pearson", "kendall", "spearman"
#'
#' @importFrom psych cortest.normal cortest.jennrich cortest.mat 
#' @importFrom stats wilcox.test kruskal.test median
#' @importFrom heplots boxM
#' @export
corrMatrix.test = function( Y, group, method = c("Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob"), method2 = c("pearson", "kendall", "spearman") ){ 

	method = match.arg( method )
	method2 = match.arg( method2 )

	if( method %in% c("Box", "Box.permute") ){

		# scale each variable, so that covariance equals correlation
		# boxM() test difference in *covariance*
		# scaling makes it test difference in *correlation*
		# Y_scale = Y
		# for( key in levels(group)){
		# 	Y_scale[group==key,] = scale(Y[group==key,])
		# }

		if( method == "Box"){
			fit = boxM( Y, group, method=method2 )
		}else{			
			fit = boxM_permute( Y, group, method=method2 )
		}

		res = list(pVal=fit$p.value, stats=fit$stat_logdet, count=NA, sign=NA)
		res$sign = ifelse( res$stats > 0, "pos", "neg")
		
		# if( nlevels(group) == 2){
		# 	# get correlations
		#   	C1 = cor( Y[group==levels(group)[1],], method=method2 )
		#   	C2 = cor( Y[group==levels(group)[2],], method=method2 )
		#   	res$stats = median(C1 - C2)
		#   	res$sign = ifelse( res$stats > 0, "pos", "neg")
		# }

	}else{

		if( nlevels(group) > 2){
			out = paste0("Method '", method, "' is only able to compare 2 groups.\nThere are ", nlevels(group), " levels in the group variable: ", paste(levels(group), collapse=', '))
			stop( out )
		}

		Y1 = Y[group==levels(group)[1],]
		Y2 = Y[group==levels(group)[2],] 

		# get correlations
	  	C1 = cor( Y1, method=method2 )
	  	C2 = cor( Y2, method=method2 )

	  	N1 = sum(group==levels(group)[1])
	  	N2 = sum(group==levels(group)[2])
		
		fit = switch( method, 
			"Steiger.fisher"= cortest.normal(C1, C2, n1=N1, n2=N2, fisher=TRUE),
			"Steiger" 		= cortest.normal(C1, C2, n1=N1, n2=N2, fisher=FALSE),
			"Jennrich"		= cortest.jennrich(C1, C2, n1=N1, n2=N2),
			"Factor" 		= cortest.mat(C1, C2, n1=N1, n2=N2),
			"Mann.Whitney"	= list(prob = 
				wilcox.test( C1[lower.tri(C1)], C2[lower.tri(C2)], paired=TRUE)$p.value),
			"Kruskal.Wallis"	= list(prob = 
				kruskal.test( C1[lower.tri(C1)], C2[lower.tri(C2)])$p.value),
			"Cai.max" = list(prob = Cai.max.test( Y1, Y2 )$pVal ),
			"Chang.maxBoot" = list(prob = Chang.maxBoot.test( Y1, Y2 )$pVal ),
			"LC.U.test" = list(prob = LC.U.test( Y1, Y2 )$pVal ),
			"WL.randProj" = list(prob = WL.randProj.test( Y1, Y2 )$pVal ),
			"Schott.Frob" = list(prob = Schott.Frob.test( Y1, Y2 )$pVal )
			)

		res = list(pVal=fit$prob, stats=median(C1 - C2), count=NA)
		res$sign = ifelse( res$stats > 0, "pos", "neg")
	}

	res
}


#' Test difference in correlation using closed form tests
#'
#' Test difference in correlation using closed form tests
#'
#' @param itObj iterator
#' @param method Specify test: "Box", "Box.permute", Steiger.fisher", "Steiger", "Jennrich", "Factor" "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob"
#' @param method2 Specify type of correlation: "pearson", "kendall", "spearman"
#'
runFastStat = function( itObj, method = c("Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis", "Cai.max", "Chang.maxBoot", "LC.U.test", "WL.randProj", "Schott.Frob"), method2 = c("pearson", "kendall", "spearman") ){ 

    method <- match.arg(method)
    method2 <- match.arg(method2)

	if( ncol(itObj$Y) >= 2 ){
	  	# test correlations
	  	res = corrMatrix.test( itObj$Y, itObj$group, method, method2)
	}else{
	  	res = list(pVal=NA, stats=NA, count=NA, sign=NA)
	}
	res$id = itObj$ID
	res$chrom = itObj$CHROM
	res$cluster = itObj$CLST
	res
}




#' Internal .evalDiffCorr
#' 
#' @param epiSignal matrix or EList of epigentic signal.  Rows are features and columns are samples
#' @param testVariable factor indicating two subsets of the samples to compare
#' @param gRanges GenomciRanges corresponding to the rows of epiSignal
#' @param clustList list of cluster assignments
#' @param npermute array of two entries with min and max number of permutations
#' @param adj.beta parameter for sLED
#' @param rho a large positive constant such that A(X)-A(Y)+diag(rep(rho,p)) is positive definite. Where p is the number of features
#' @param sumabs.seq sparsity parameter
#' @param BPPARAM parameters for parallel evaluation
#' @param method "sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney" or "Kruskal.Wallis" 
#' @param method2 Specify type of correlation: "pearson", "kendall", "spearman"
#' 
#' 
#' @return list of result by chromosome and clustList
#'
#' @import BiocParallel
#' @import foreach
#' @importFrom progress progress_bar
#' @importFrom data.table data.table
.evalDiffCorr = function(epiSignal, testVariable, gRanges, clustList, npermute = c(100, 10000, 1e5), adj.beta=0, rho=0, sumabs.seq=1, BPPARAM = SerialParam(), method=c("sLED", "Box", "Box.permute", "Steiger.fisher", "Steiger", "Jennrich", "Factor", "Mann.Whitney", "Kruskal.Wallis"), method2=c("pearson", "kendall", "spearman")){

	if( nrow(epiSignal) != length(gRanges)){
		stop("Number of rows in epiSignal must equal number of entries in gRanges")
	}
	if( length(testVariable) != ncol(epiSignal) ){
		stop("Number of columns in epiSignal must equal number of entries in testVariable")
	}
	if( ! is(testVariable, 'factor') ){
		stop("Entries in testVariable must be a factor with entries in two levels")
	}
	method = match.arg( method )
	if( method == "sLED" && (length(npermute) < 2 || is.unsorted(npermute) ) ){
		stop("npermute must have 2 or 3 increasing entries")
	}
	testVariable = droplevels(testVariable)
	if( (nlevels(testVariable) < 2) || (nlevels(testVariable) > 2 && ! (method %in% c("Box", "Box.permute"))) ){
		out = paste0("Method '", method, "' is only able to compare 2 groups.\nThere are ", nlevels(testVariable), " levels in the group variable: ", paste(levels(testVariable), collapse=', '))
		stop( out )
	}
	if( min(table(testVariable)) < 10){
		stop("Need at last 10 samples in smallest class of testVariable")
	}	
	if( any(!unique(unlist(lapply(clustList, names))) %in% seqnames(gRanges)@values) ){
		stop("Chromosomes from clustList must be in gRanges")
	} 
	if( length(adj.beta) > 1 ){
		stop("adj.beta must be a scalar")
	}
	if( length(sumabs.seq) > 1 ){
		stop("sumabs.seq must be a scalar")
	}

	method2 = match.arg( method2 )

	if( method2 != "pearson" && method == "sLED"){
		stop("sLED only evalutes pearson correlation")
	}

	# allClusters = unlist(lapply(names(clustList), function(x) paste0(x, '_', unique(clustList[[x]]))))

	# Divide data into two sets
	set1 = which(testVariable == levels(testVariable)[1])
	set2 = which(testVariable == levels(testVariable)[2])

	# create single data.table with cluster and peak info
	dfClust = lapply(names(clustList), function(id){
		res = lapply( names(clustList[[id]]), function(chrom){
			data.frame(id=id, chrom=chrom, cluster=clustList[[id]][[chrom]], peak=names(clustList[[id]][[chrom]]), stringsAsFactors=FALSE)
		})
		do.call("rbind", res)
	})
	dfClust = data.table(do.call('rbind', dfClust))

	# check size of clusters
	cat("Note that clusters of 2 or fewer features are omitted from analysis\n\n")

	# only get peaks that are in the the epiSignal dataset
	peak = NA
	dfClust = dfClust[peak %in% rownames(epiSignal),]

	# get unique chrom and size of cluster
	.N = NA
	dfClustCounts = dfClust[,.N, by=c("id", "chrom", "cluster")]

	# sort by size within chromosomes
	# evaluate sLED starting with largest gene set decreasing
	# 	this makes time estimate more accurate
	# dfClustCountsSort = dfClustCounts[,.SD[order(N, decreasing=TRUE),],by="chrom"]
	# sort by size across chromosomes
	.SD = N = NA
	dfClustCountsSort = dfClustCounts[,.SD[order(N, decreasing=TRUE),]]

	n_clusters = nrow(dfClustCountsSort)

	cat("# Clusters:", n_clusters, '\n')

	# Evaluate statistics with permutations
	#######################################

	# res = mclapply( seq_len(nrow(dfClustUnique)), runSled, dfClustUnique, dfClust, epiSignal, set1, set2, npermute=c(10, 100), mc.cores=1)

	# lapply(1, runSled, dfClustUnique, dfClust, epiSignal, set1, set2, npermute)

	 # runSled(1, dfClustUnique, dfClust, epiSignal, set1, set2, npermute)

	# combinedResults = bplapply( seq_len(nrow(dfClustUnique)), runSled, dfClustUnique, dfClust, epiSignal, set1, set2, npermute, BPPARAM=BPPARAM)
	
	if( method == "sLED"){
		cat("Initial pass through all clusters...\n")
		# run with iterators
		it = clustIter( dfClustCountsSort, dfClust, epiSignal, testVariable )
		
		combinedResults = bpiterate( it$nextElem, runSled2, npermute, adj.beta, rho, sumabs.seq, SerialParam(), BPPARAM=BPPARAM)

		# create data.frame
		df = list()
		df$permCounts = vapply(combinedResults, function(x) x$count, numeric(1))
		df$id = vapply(combinedResults, function(x) x$id, 'character')
		df$chromArray = vapply(combinedResults, function(x) x$chrom, 'character')
		df$clustArray = vapply(combinedResults, function(x) x$cluster, numeric(1))
		df = data.table(data.frame(df, stringsAsFactors=FALSE))

		numPassCutoff = sum(df$permCounts < 10, na.rm=TRUE)

		if( numPassCutoff > 0 ){
			cat("Intensive second pass...\n")

			# parallelize run of each cluster
			pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",	total = numPassCutoff, width= 60, clear=FALSE)

			itGlobal = clustIter( dfClustCountsSort, dfClust, epiSignal, testVariable )
			
			npermute2 = sort(c(npermute[2]*10, npermute[3]))

			count = 0
			while( ! is.null( it <- itGlobal$nextElem() ) ){ 

				# fit location in df of the features set in it
				i = which( (df$id == it$ID) & (df$chromArray == it$CHROM) & (df$clustArray == it$CLST))

				prmCnt = df$permCounts[i]
				# if permCounts is too small, run intensive parallel analysis
				if( !is.na(prmCnt) & (prmCnt< 10) & (length(npermute) >=3) ){

					# update progress bar
					pb$update( min(count / numPassCutoff, 1) )
					count = count + 1

					# run analysis
					# suppress progressbar
					suppressMessages({
					combinedResults[[i]] <- runSled2( it,npermute2, adj.beta, rho, sumabs.seq, BPPARAM )
					})
				}
			}
			pb$update( 1.0 )
		}
	}else{

		# run with iterators
		it = clustIter( dfClustCountsSort, dfClust, epiSignal, testVariable )
		
		combinedResults = bpiterate( it$nextElem, runFastStat, method, method2, BPPARAM=BPPARAM)
	}

	# return list of lists
	######################

	cat("Combining results...\n")

	resList = list()
	chromArray = unique(dfClustCounts$chrom)
	idArray = unique(dfClustCounts$id)

	# unique(vapply(combinedResults, function(x) x$chrom, "character"))
	for( chrom in chromArray){

		key = vapply(combinedResults, function(x) paste(x$id, x$cluster), 'character')

		idx = vapply(combinedResults, function(x){ 
			(paste(x$id, x$cluster) %in% key) && (x$chrom == chrom)
			}, logical(1)) 

		resList[[chrom]] = combinedResults[idx]
		names(resList[[chrom]]) = vapply(resList[[chrom]], function(x) paste(x$id, x$cluster), 'character')
	}
	names(resList) = chromArray
	
	cat("\n")
	new("sLEDresults", resList)
}





#' @import sLED
.sLED = function (X, Y, adj.beta=0, rho = 0, sumabs.seq = 1,
    npermute = 100, seeds = NULL,
    verbose = TRUE, niter = 20, trace = FALSE, BPPARAM=SerialParam()){
    D.hat <- sLED:::getDiffMatrix(X, Y, adj.beta)

    pma.results <- sLED:::sLEDTestStat(Dmat = D.hat, rho = rho, sumabs.seq = sumabs.seq, niter = niter, trace = trace)
    Tn <- pma.results$stats
    n1 <- nrow(X)
    n2 <- nrow(Y)
    Z <- rbind(X, Y)
    permute.results <- .sLEDpermute(Z = Z, n1 = n1, n2 = n2, adj.beta = adj.beta,
        rho = rho, sumabs.seq = sumabs.seq, npermute = npermute,
        seeds = seeds, verbose = verbose,
        niter = niter, trace = trace, BPPARAM=BPPARAM)

    count = rowSums(permute.results$Tn.permute > Tn)
    pVal = (count+1) / (npermute+1)
    return(c(pma.results, permute.results, list(Tn = Tn, count=count, pVal = pVal)))
}


#' @import sLED
#' @importFrom BiocParallel bplapply
.sLEDpermute = function (Z, n1, n2, adj.beta=0, rho = 0, sumabs.seq = 1,
    npermute = 100, seeds = NULL,
    verbose = TRUE, niter = 20, trace = FALSE, BPPARAM=SerialParam()){
    if (verbose) {
        cat(npermute, "permutation started:\n")
    }
            
    perm.results <- bplapply(seq_len(npermute), sLED:::sLEDOnePermute,
        Z = Z, n1 = n1, n2 = n2, seeds = seeds, sumabs.seq = sumabs.seq,
        adj.beta = adj.beta, rho = rho, verbose = FALSE,
        niter = niter, trace = trace, BPPARAM=BPPARAM)   
    
    ntest <- length(sumabs.seq)
    Tn.permute <- matrix(NA, nrow = ntest, ncol = npermute)
    Tn.permute.sign <- matrix(NA, nrow = ntest, ncol = npermute)
    for (i in seq_len(npermute)) {
        Tn.permute[, i] <- perm.results[[i]]$Tn.permute
        Tn.permute.sign[, i] <- perm.results[[i]]$Tn.permute.sign
    }
    if (verbose) {
        cat("permutations finished.", fill = TRUE)
    }
    return(list(Tn.permute = Tn.permute, Tn.permute.sign = Tn.permute.sign))
}










