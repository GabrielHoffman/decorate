
test_boxM = function(){

		
	Y = as.matrix(iris[,1:4])
	group = iris[,5]

	method = "pearson"
	res1 = boxM( Y, group, method=method)
	res2 = decorate:::boxM_fast( Y, group, method=method)

	method = "pearson"
	res3 = boxM( Y, group, method=method)
	res4 = decorate:::boxM_fast( Y, group, method=method)

	checkEquals( as.numeric(res1$statistic), res2$X2) && 
	checkEquals( as.numeric(res3$statistic), res4$X2)
}
