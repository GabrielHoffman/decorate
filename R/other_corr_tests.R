
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
#' 
#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#' 
#' @return list of p-value, estimate and method used 
#' 
#' @details The statistical test used depends on the variable specified.
#' if variable is factor with multiple levels, use Kruskal-Wallis test
#' if variable is factor with 2 levels, use Wilcoxon test
#' if variable is continuous, use Wilcoxon test
#' 
#' @examples
#' # load iris data
#' data(iris)
#' 
#' # variable is factor with multiple levels
#' # use kruskal.test
#' delaneau.test( iris[,1:4], iris[,5] )
#' 
#' # variable is factor with 2 levels
#' # use wilcox.test
#' delaneau.test( iris[1:100,1:4], iris[1:100,5] )
#' 
#' # variable is continuous
#' # use cor.test with spearman
#' delaneau.test( iris[,1:4], iris[,1] )
#' 
#' @importFrom stats cor.test
#' @export
#' @seealso delaneau.score sle.test
delaneau.test = function( Y, variable, method = c("pearson", "kendall", "spearman") ){

      # compute sample level scores
      score = delaneau.score( Y, method )

       # Test of association between score and variable
      if( is.factor(variable) ){
            variable = droplevels(variable)
            if( nlevels(variable) == 2 ){
                  x = score[variable==levels(variable)[1]]
                  y = score[variable==levels(variable)[2]]
                  fit = wilcox.test( x, y ) 
                  result = list(p.value = fit$p.value, estimate = median(y) - median(x), method="wilcox.test")  
            }else{
                  fit = kruskal.test( score, variable )
                  result = list(p.value = fit$p.value, estimate = fit$statistic, method="kruskal.test")   
            }
      }else{
            fit = cor.test( score, variable, method="spearman") 
            result = list(p.value = fit$p.value, estimate = fit$estimate, method="cor.test") 
      }

      result
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
sle.score = function( Y, method = c("pearson", "kendall", "spearman"), rho=.05, sumabs=1 ){
      
      method = match.arg( method )

      # correlation for all samples
      C_all = cor(Y, method=method)

      # Compute leave-one-out correlations and statistics
      sampleScore = vapply( seq_len(nrow(Y)), function( i ){

            # correlation dropping sample i
            C_i = cor( Y[-i,], method=method)

            # Difference matrix
            D.mat = C_all - C_i

            # if more thatn 2 features
            if( nrow(D.mat) > 2){
                  # Evaluate sLED statistic
                  res = sLED:::sLEDTestStat( D.mat, rho=rho, sumabs.seq = sumabs)
            }else{
                  # Evaluate statistic without sparsity
                  p = nrow(D.mat)
                 
                  d.pos = eigen(D.mat + rho*diag(p))$values[1]
                  d.neg = eigen(-1*D.mat + rho*diag(p))$values[1]

                  res = list(stats = ifelse( d.pos > d.neg, d.pos - rho, -1*(d.neg - rho)))
            }

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
#' Score impact of each sample on sparse leading eigen-value and then peform test of association with variable using non-parametric test
#'
#' @param Y data matrix with samples on rows and variables on columns
#' @param variable variable with number of entries must equal nrow(Y).  Can be discrete or continuous.
#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#' @param rho a positive constant such that cor(Y) + diag(rep(rho,p)) is positive definite.
#' @param sumabs regularization paramter. Value of 1 gives no regularization, sumabs*sqrt(p) is the upperbound of the L_1 norm of v,controling the sparsity of solution. Must be between 1/sqrt(p) and 1.
#' 
#' @return list of p-value, estimate and method used 
#' 
#' @details The statistical test used depends on the variable specified.
#' if variable is factor with multiple levels, use Kruskal-Wallis test
#' if variable is factor with 2 levels, use Wilcoxon test
#' if variable is continuous, use Wilcoxon test
#' 
#' @examples
#' # load iris data
#' data(iris)
#' 
#' # variable is factor with multiple levels
#' # use kruskal.test
#' sle.test( iris[,1:4], iris[,5] )
#' 
#' # variable is factor with 2 levels
#' # use wilcox.test
#' sle.test( iris[1:100,1:4], iris[1:100,5] )
#' 
#' # variable is continuous
#' # use cor.test with spearman
#' sle.test( iris[,1:4], iris[,1] )
#' 
#' @importFrom stats cor.test
#' @export
#' @seealso sle.score delaneau.test
sle.test = function( Y, variable, method = c("pearson", "kendall", "spearman"), rho=0, sumabs=1 ){

      # compute sample level scores
      score = sle.score( Y, method )

       # Test of association between score and variable
      if( is.factor(variable) ){
            variable = droplevels(variable)
            if( nlevels(variable) == 2 ){
                  x = score[variable==levels(variable)[1]]
                  y = score[variable==levels(variable)[2]]
                  fit = wilcox.test( x, y ) 
                  result = list(p.value = fit$p.value, estimate = median(y) - median(x), method="wilcox.test")  
            }else{
                  fit = kruskal.test( score, variable )
                  result = list(p.value = fit$p.value, estimate = fit$statistic, method="kruskal.test")   
            }
      }else{
            fit = cor.test( score, variable, method="spearman") 
            result = list(p.value = fit$p.value, estimate = fit$estimate, method="cor.test") 
      }

      result
}












