
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

			# get two subsets of datas
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


clustIterBatch = function( dfClust, epiSignal, testVariable, n_chunks = 100 ){

	n_clusters = countClusters(dfClust)

	# divide clusters into batches
	idx <- parallel::splitIndices(n_clusters, min(n_clusters, n_chunks))

    i <- 0L
    f = function() {
        if (i == length(idx)){
            return(NULL)
        }
        i <<- i + 1L

        # get subset of clusters
        clustObj = getClusterSubset( dfClust, idx[[i]] )

        # get cluster ids of the retained clusters
        df = attr(clustObj, "getClusterNames")
       
       	# extract peak ids
        peakIds = lapply( seq_len(nrow(df)), function(i){
        	getFeaturesInCluster( clustObj, df$chrom[i], as.character(df$cluster[i]), as.character(df$parameter[i]))
        	})
        peakIds = unlist(peakIds)

        list( 	clustObj = clustObj,
        		signalObj = epiSignal[peakIds,,drop=FALSE],
        		testVariable = testVariable )
    }

    attr( f, "n_chunks") = length(idx)
    f
}


corrIterBatch = function( treeList, treeListClusters, df_count, n_chunks = 100 ){

    index = meanClusterSize = NA

    df_count = data.table(df_count)
    df_count[,index:=seq_len(nrow(df_count))]

    # divide batches to include a single chrom and meanClusterSize
    idx = list()
    i = 1
    for(chrom_ in unique(df_count$chrom)){
        for(id in unique(df_count$meanClusterSize)){
            idx[[i]] = df_count[chrom==chrom_ & meanClusterSize==id, index]
            i = i + 1
        }
    }

    i <- 0L
    f = function() {
        if (i == length(idx)){
            return(NULL)
        }
        i <<- i + 1L

        df_sub = df_count[idx[[i]],]

        chrom = df_sub[,unique(chrom)]
        id = as.character(df_sub[,unique(meanClusterSize)])

        treeList_sub = treeList[[chrom]]
        treeListClusters_sub = treeListClusters[[id]][[chrom]]

        list(   treeList            = treeList_sub,
                treeListClusters    = treeListClusters_sub,
                df_count            = df_sub, 
                chrom               = chrom,
                id                   = id)         
    }

    attr( f, "n_chunks") = length(idx)
    f
}












