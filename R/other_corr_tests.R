
#' Score impact of each sample on correlation sturucture 
#'
#' Score impact of each sample on correlation sturucture.  Compute correlation using all samples (i.e. C), then compute correlation omitting sample i (i.e. Ci).  The score the sample i is based on the difference between C and Ci.
#'
#' @param Y data matrix with samples on rows and variables on columns
#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#' 
#' @return score for each sample measure impact on correlation structure 
#' 
#' @examples
#' # load iris data
#' data(iris)
#' 
#' # Evalaute score on each sample
#' delaneau.score( iris[,1:4] )
#' 
#' @importFrom stats cor
#' @export
#' @seealso delaneau.test
delaneau.score = function( Y, method = c("pearson", "kendall", "spearman") ){

      method = match.arg( method )

      # correlation for all samples
      C_all = cor(Y, method=method)

      # mean correlation
      mu_c = mean(C_all[lower.tri(C_all)]) 

      # Compute leave-one-out correlations and statistics
      sampleScore = vapply( seq_len(nrow(Y)), function( i ){
            C_i = cor( Y[-i,], method=method)

            mean(C_i[lower.tri(C_i)]) - mu_c
            }, numeric(1))

      if( !is.null(rownames(Y)) ){
            names(sampleScore) = rownames(Y)
      }
     sampleScore      
}


#' Test association between correlation sturucture and variable
#'
#' Score impact of each sample on correlation sturucture and then peform test of association with variable using Kruskal-Wallis test
#'
#' @param Y data matrix with samples on rows and variables on columns
#' @param variable variable with number of entries must equal nrow(Y).  Can be discrete or continuous.

#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#' 
#' @return score for each sample measure impact on correlation structure 
#' 
#' @examples
#' # load iris data
#' data(iris)
#' 
#' # Evalaute score on each sample
#' delaneau.test( iris[,1:4], iris[,5] )
#' 
#' @importFrom stats cor
#' @export
#' @seealso delaneau.score
delaneau.test = function( Y, variable, method = c("pearson", "kendall", "spearman") ){

      # compute sample level scores
      score = delaneau.score( Y, method )

      # Test of association between score and variable
      kruskal.test( score, variable )   
}





#' Score impact of each sample on sparse leading eigen-value
#'
#' Score impact of each sample on sparse leading eigen-value.  Compute correlation using all samples (i.e. C), then compute correlation omitting sample i (i.e. Ci).  The score the sample i is based on sparse leading eigen-value of the diffrence between C and Ci.
#'
#' @param Y data matrix with samples on rows and variables on columns
#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#' @param rho a positive constant such that cor(Y) + diag(rep(rho,p)) is positive definite.
#' @param sumabs regularization paramter. Value of 1 gives no regularization, sumabs*sqrt(p) is the upperbound of the L_1 norm of v,controling the sparsity of solution. Must be between 1/sqrt(p) and 1.
#' 
#' @return score for each sample measure impact on correlation structure 
#' 
#' @examples
#' # load iris data
#' data(iris)
#' 
#' # Evalaute score on each sample
#' sle.score( iris[,1:4] )
#' 
#' @importFrom stats cor
#' @export
#' @seealso sle.test
sle.score = function( Y, method = c("pearson", "kendall", "spearman"), rho=.1, sumabs=1 ){
      method = match.arg( method )

      # correlation for all samples
      C_all = cor(Y, method=method)

      # Compute leave-one-out correlations and statistics
      sampleScore = vapply( seq_len(nrow(Y)), function( i ){

            # correlation dropping sample i
            C_i = cor( Y[-i,], method=method)

            # Difference matrix
            D.mat = C_all - C_i

            # Evaluate sLED statistic
            res = sLED:::sLEDTestStat( D.mat, rho=rho, sumabs.seq = sumabs)

            # return test stat
            res$stats
            }, numeric(1))

      if( !is.null(rownames(Y)) ){
            names(sampleScore) = rownames(Y)
      }
      sampleScore   
}


#' Test association between sparse leading eigen-value and variable
#'
#' Score impact of each sample on sparse leading eigen-value and then peform test of association with variable using Kruskal-Wallis test
#'
#' @param Y data matrix with samples on rows and variables on columns
#' @param variable variable with number of entries must equal nrow(Y).  Can be discrete or continuous.
#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#' @param rho a positive constant such that cor(Y) + diag(rep(rho,p)) is positive definite.
#' @param sumabs regularization paramter. Value of 1 gives no regularization, sumabs*sqrt(p) is the upperbound of the L_1 norm of v,controling the sparsity of solution. Must be between 1/sqrt(p) and 1.
#' 
#' @return score for each sample measure impact on correlation structure 
#' 
#' @examples
#' # load iris data
#' data(iris)
#' 
#' # Evalaute score on each sample
#' sle.test( iris[,1:4], iris[,5] )
#' 
#' @importFrom stats cor
#' @export
#' @seealso sle.score
sle.test = function( Y, variable, method = c("pearson", "kendall", "spearman"), rho=0, sumabs=1 ){

       # compute sample level scores
      score = sle.score( Y, method )

      # Test of association between score and variable
      kruskal.test( score, variable )   
}



















