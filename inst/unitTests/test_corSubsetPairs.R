
#' @export RUnit
test_corSubsetPairs = function(){

	N = 1000
	Y = matrix(rnorm(N*100), 100, N)

	v1 = sample.int(N, 2000, replace=TRUE)
	v2 = sample.int(N, 2000, replace=TRUE)

	M = corSubsetPairs(t(Y), v1,v2)

	# compute correlation using corSubsetPairs
	cor1 = M[v1[1], v2[1]]

	# direct compute correlation
	cor2 = cor(Y[,v1[1]], Y[,v2[1]])

	checkEquals(cor1, cor2)
}
