#' Compute correlations between pairs of features
#'
#' Compute correlations between pairs of features given in idxi and idxj
#'
#' @param Y matrix where rows are features
#' @param idxi indecies
#' @param idxj indecies
#' @param method specify which correlation method: "pearson" or "spearman"
#' @param silent suppress messages 
#' @param setNANtoZero replace NAN correlation values with a zero
#'
#' @return 
#' Compute local correlations between for all k: cor(Y[,idxi[k]], Y[,idxj[k]])
#' 
#' @examples
#' # Simulate simple dataset
#' N = 600
#' Y = matrix(rnorm(N*100), 100, N)
#' 
#' # select pairs to compute correlations between
#' i1 = sample.int(N, 200, replace=TRUE)
#' i2 = sample.int(N, 200, replace=TRUE)
#' 
#' # evaluate all piars
#' C = corSubsetPairs(t(Y), i1,i2)
#' 
#' # show value
#' C[i1[10], i2[10]]
#' 
#' # show values from evaluating this pair directly
#' cor(Y[,i1[10]], Y[,i2[10]])
#' 
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @importFrom Matrix sparseMatrix
#' @importFrom utils object.size
#' @export
corSubsetPairs <- function(Y, idxi, idxj, method = c("pearson", "spearman"), silent=FALSE, setNANtoZero=FALSE) {

	# check function arguments

	# lengths of indecies must be equal
	if( length(idxi) != length(idxj) ){
		stop("Lengths of idxi and idxj must be equal: ", length(idxi) ,' != ', length(idxj))
	}

	# max index cannot be larger than ncol(Y)
	idx_range = range(c(idxi, idxj))
	if( idx_range[2] > nrow(Y) ){
		stop("Entry indxi or idxj exceeeds nrow(Y): ", idx_range[2] ,' > ', nrow(Y))
	}

	# min index must be greater than zero
	if( idx_range[1] <= 0 ){
		stop("Entry in indxi or idxj is <= 0: ", idx_range[1])
	}

	method = match.arg( method )

	# Result is symmatric upper triangle, 
	# so any entries where i > j should be flipped
	# all pairs should be unique
	# df_idx = unique(t(apply(cbind(idxi, idxj), 1, function(x) sort(x))))
	# ============
	# sort
	vi = pmin(idxi, idxj)
	vj = pmax(idxi, idxj)

	# unique
	# *much* faster than using cbind and comparing rows of matrix
	unq = which(duplicated(paste(vi, vj)))
	if( length(unq) > 0){
		vi = vi[-unq]
		vj = vj[-unq]
	}

	if( method == "spearman"){
		Rank <- function(u) {
	        if (length(u) == 0L)
	            u
	        else if (is.matrix(u)) {
	            if (nrow(u) > 1L)
	                apply(u, 2L, rank, na.last = "keep")
	            else row(u)
	        }
	        else rank(u, na.last = "keep")
	    }

		Y = scale(Rank( t(Y) ))
	}else{
		Y = scale(t(Y))
	}


	# only computes crossprod(y1, y2)
	# This is only the correlation value if mean of cols is zero and 
	# sd == 1
    rho = as.numeric(.Call('_decorate_corSubsetPairs', PACKAGE = 'decorate', Y, vi, vj))

    N = ncol(Y)
	sparsity = length(idxi) / N^2

	# Sparse symmetric matrix
	#   create symmetrics covariance matrix
	#   using only upper triangle
	#=======================================

	if( setNANtoZero ){
		idx = is.nan(rho)
		if( any(idx) ){
			rho[which(idx)] = 0
		}
	}

	# Set diagonals to 1
	vi = append(vi, seq_len(N))
	vj = append(vj, seq_len(N))
	rho = append(rho, rep(1, N))

	# if use.last.ij == FALSE, values at repeated indeces are summed
	# if TRUE, only use last value
	M = sparseMatrix( i = vi,
                j = vj,
                x = rho,
                dimnames = list(colnames(Y), colnames(Y)),
                dims = c(N,N), symmetric=TRUE)

	if( ! silent ){
		message("Covariance matrix...\n")
		message(" sparsity:", sprintf("%.3f%s", 100 - sparsity*100, ' %\n'))
		message(" memory usage:", format(object.size(M), "MB"), '\n')
	}

	M
}


#' Box's M-test
#'
#' boxM performs the Box's (1949) M-test for homogeneity of covariance matrices obtained from multivariate normal data according to one or more classification factors. The test compares  the product of the log determinants of the separate covariance  matrices to the log determinant of the pooled covariance matrix,    analogous to a likelihood ratio test. The test statistic uses a chi-square approximation. Uses permutations to estimate the degrees of freedom under the null
#' @param Y response variable matrix
#' @param group a factor defining groups, number of entries must equal nrow(Y)
#' @param method Specify type of correlation: "pearson", "spearman"
#' 
#' 
#' @examples
#' data(iris)
#'
#' boxM_fast( as.matrix(iris[, 1:4]), iris[, "Species"])
#'
#' @seealso heplots::boxM
#' @export
boxM_fast = function( Y, group, method=c("pearson", "spearman")){

	method = match.arg( method )

	out = .Call('_decorate_boxM_fast', PACKAGE = 'decorate', Y, group, method)

	out$logdet = as.numeric(out$Si_logDet)
	out$stat_logdet = ifelse(length(out$logdet) == 2, out$logdet[2] - out$logdet[1], NA)
	out
}



