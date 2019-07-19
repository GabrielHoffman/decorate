
#' @export RUnit
test_corSubsetPairs = function(){

	N = 1000
	Y = matrix(rnorm(N*100), 100, N)

	v1 = sample.int(N, 2000, replace=TRUE)
	v2 = sample.int(N, 2000, replace=TRUE)

	method = "pearson"
	M = corSubsetPairs(t(Y), v1,v2, method)
	# compute correlation using corSubsetPairs
	cor1 = M[v1[1], v2[1]]
	# direct compute correlation
	cor2 = cor(Y[,v1[1]], Y[,v2[1]], method=method)

	method = "spearman"
	M = corSubsetPairs(t(Y), v1,v2, method)
	# compute correlation using corSubsetPairs
	cor3 = M[v1[1], v2[1]]
	# direct compute correlation
	cor4 = cor(Y[,v1[1]], Y[,v2[1]], method=method)


	checkEquals(cor1, cor2) & checkEquals(cor3, cor4)
}
