
# Compute log determinate based on only the first N eigen values
det_partial = function( C, tol=1e-10 ){

   if( ncol(C) != nrow(C) ){
      stop("Matrix must be square")
   }

   dcmp = eigen( C, only.values=TRUE)

   evalues = dcmp$values

   # log determinant is sum of log eigen-values
   # partial sum of 1:N values
   sum( log(evalues[evalues > tol]) )
}

#' Box's M-test 
#' 
#' Box's M-test 
#' 
#' @param Y response matrix
#' @param group factor defining groups
#' @param tol tolerance for eigen values
#' @param fxn define function.  Here default is cor to compare correlation structure.  Use cov to compare covariance structure like in  heplots::boxM  
#' @param method specify which correlation method: "pearson", "kendall" or "spearman"
#'
#' @examples
#' data(iris)
#' 
#' boxM( iris[,1:4], iris[,5])
#' 
#' @importFrom stats complete.cases aggregate


#' @export
boxM <- function(Y, group, tol=1e-10, fxn=cor, method= c("pearson", "kendall", "spearman")){
   dname <- deparse(substitute(Y))
   method = match.arg( method )

   if (!inherits(Y, c("data.frame", "matrix")))
      stop(paste(dname, "must be a numeric data.frame or matrix!"))
   if (length(group) != nrow(Y))
      stop("incompatible dimensions in Y and group!")

   Y <- as.matrix(Y)
   gname <- deparse(substitute(group))
   if (!is.factor(group)) group <- as.factor(as.character(group))

   valid <- complete.cases(Y, group)
   if (nrow(Y) > sum(valid)){ 
      warning(paste(nrow(Y) - sum(valid)), " cases with missing data have been removed.")
   }
   Y <- Y[valid,]
   group <- group[valid]

   p <- ncol(Y)
   nlev <- nlevels(group)
   lev <- levels(group)
   dfs <- tapply(group, group, length) - 1
   if( any(dfs < p) ){ 
      warning("there are one or more levels with less observations than variables!")
   }

   mats <- aux <- list()
   for(i in seq_len(nlev) ){
      mats[[i]] <-  cor( Y[group == lev[i], ], method=method )  
      aux[[i]] <- mats[[i]] * dfs[i]
   }
   names(mats) <- lev

   pooled <- Reduce("+", aux) / sum(dfs)

   # rank = min( rank, ncol(Y) )

   # logdet <- log(unlist(lapply(mats, det)))
   logdet <- vapply(mats, det_partial, tol=tol, numeric(1))

   # logdet_pooled = log(det(pooled))
   logdet_pooled = det_partial( pooled, tol=tol )

   minus2logM <- sum(dfs) * logdet_pooled - sum(logdet * dfs)
   sum1 <- sum(1 / dfs) 

   Co <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
     (nlev - 1))) * (sum1 - (1 / sum(dfs)))

   X2 <- minus2logM * (1 - Co)

   dfchi <- (choose(p, 2) + p) * (nlev - 1)

   pval <- pchisq(X2, dfchi, lower.tail = FALSE)

   means <- aggregate(Y, list(group), mean)
   rn <- as.character(means[,1])
   means <- means[,-1]
   means <- rbind( means, colMeans(Y) )
   rownames(means) <- c(rn, "pooled")

   logdet <- c(logdet, pooled=log(det(pooled)))
   df <- c(dfs, pooled=sum(dfs))
   out <- structure(
      list(statistic = c("Chi-Sq (approx.)" = X2),
         parameter = c(df = dfchi),
         p.value = pval,
         cov = mats, pooled = pooled, logDet = logdet, means = means, df=df,
         data.name = dname, group = gname,
         method = "Box's M-test for Homogeneity of Covariance Matrices"
         ),
      class = c("htest", "boxM")
      )
   out$logdet = logdet
   out$stat_logdet = ifelse(length(logdet) == 3, logdet[2] - logdet[1], NA)
   return(out)
}


# Compute boxM result using only the first *rank* eigen-values 
boxM_basic <- function(Y, group, tol=1e-10, fxn=cor, method= c("pearson", "kendall", "spearman")){
   # dname <- deparse(substitute(Y))
   # method = match.arg( method )

   # if (!inherits(Y, c("data.frame", "matrix")))
   #    stop(paste(dname, "must be a numeric data.frame or matrix!"))
   # if (length(group) != nrow(Y))
   #    stop("incompatible dimensions in Y and group!")

   # Y <- as.matrix(Y)
   # gname <- deparse(substitute(group))
   # if (!is.factor(group)) group <- as.factor(as.character(group))

   # valid <- complete.cases(Y, group)
   # if (nrow(Y) > sum(valid)){ 
   #    warning(paste(nrow(Y) - sum(valid)), " cases with missing data have been removed.")
   # }
   # Y <- Y[valid,]
   # group <- group[valid]

   p <- ncol(Y)
   nlev <- nlevels(group)
   lev <- levels(group)
   dfs <- tapply(group, group, length) - 1
   # if( any(dfs < p) ){ 
   #    warning("there are one or more levels with less observations than variables!")
   # }

   mats <- aux <- list()
   for(i in seq_len(nlev) ){
      mats[[i]] <-  cor( Y[group == lev[i], ], method=method )  
      aux[[i]] <- mats[[i]] * dfs[i]
   }
   names(mats) <- lev

   pooled <- Reduce("+", aux) / sum(dfs)

   # rank = min( rank, ncol(Y) )

   # logdet <- log(unlist(lapply(mats, det)))
   logdet <- vapply(mats, det_partial, tol=tol, numeric(1))

   # logdet_pooled = log(det(pooled))
   logdet_pooled = det_partial( pooled, tol=tol )

   minus2logM <- sum(dfs) * logdet_pooled - sum(logdet * dfs)
   sum1 <- sum(1 / dfs) 

   Co <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
     (nlev - 1))) * (sum1 - (1 / sum(dfs)))

   X2 <- minus2logM * (1 - Co)

   dfchi <- (choose(p, 2) + p) * (nlev - 1)

   pval <- pchisq(X2, dfchi, lower.tail = FALSE)

   # means <- aggregate(Y, list(group), mean)
   # rn <- as.character(means[,1])
   # means <- means[,-1]
   # means <- rbind( means, colMeans(Y) )
   # rownames(means) <- c(rn, "pooled")

   # logdet <- c(logdet, pooled=log(det(pooled)))
   # df <- c(dfs, pooled=sum(dfs))
   # out <- structure(
   #    list(statistic = c("Chi-Sq (approx.)" = X2),
   #       parameter = c(df = dfchi),
   #       p.value = pval,
   #       cov = mats, pooled = pooled, logDet = logdet, means = means, df=df,
   #       data.name = dname, group = gname,
   #       method = "Box's M-test for Homogeneity of Covariance Matrices"
   #       ),
   #    class = c("htest", "boxM")
   #    )
   out = list( statistic = X2, p.value = pval, parameter = c(df = dfchi))
   out$logdet = logdet
   out$stat_logdet = ifelse(length(logdet) == 3, logdet[2] - logdet[1], NA)

   return(out)
}




#' Box's M-test
#'
#' boxM performs the Box's (1949) M-test for homogeneity of covariance matrices obtained from multivariate normal data according to one or more classification factors. The test compares  the product of the log determinants of the separate covariance  matrices to the log determinant of the pooled covariance matrix,    analogous to a likelihood ratio test. The test statistic uses a chi-square approximation. Uses permutations to estimate the degrees of freedom under the null
#' @param Y response variable matrix
#' @param group a factor defining groups, number of entries must equal nrow(Y)
#' @param nperm number of permutations of group variable used to estimate degrees of freedom under the null
#' @param method Specify type of correlation: "pearson", "kendall", "spearman"
#' 
#' 
#' @examples
#' data(iris)
#'
#' boxM_permute(iris[, 1:4], iris[, "Species"])
#'
#' @return list of p.value, test statistic, and df.approx estimated by permutation
#'
#' @importFrom stats optimize pchisq dchisq
#' @seealso heplots::boxM
#' @export
boxM_permute = function(Y, group, nperm=200, method=c("pearson", "kendall", "spearman")){

   dname <- deparse(substitute(Y))
   
   if (!inherits(Y, c("data.frame", "matrix")))
      stop(paste(dname, "must be a numeric data.frame or matrix!"))
   if (length(group) != nrow(Y))
      stop("incompatible dimensions in Y and group!")

   Y <- as.matrix(Y)
   gname <- deparse(substitute(group))
   if (!is.factor(group)) group <- as.factor(as.character(group))

   valid <- complete.cases(Y, group)
   if (nrow(Y) > sum(valid)){ 
      warning(paste(nrow(Y) - sum(valid)), " cases with missing data have been removed.")
   }
   Y <- Y[valid,]
   group <- group[valid]

   p <- ncol(Y)
   nlev <- nlevels(group)
   lev <- levels(group)
   dfs <- tapply(group, group, length) - 1
   if( any(dfs < p) ){ 
      warning("there are one or more levels with fewer observations than variables!")
      return(list( p.value = NA,
      statistic = NA,
      df.approx = NA,
      stat_logdet = NA))
   }

   method = match.arg( method )  

   # fit statistic for real data
   fit = boxM_fast( Y, group, method)

   # fit statistic for permutated data
   # get chisq statistic
   chisq_stats = vapply( seq_len(nperm), function(i){
      idx = sample.int(length(group), length(group))
      fit = boxM_fast( Y, group[idx], method)
      fit$X2
      }, numeric(1))

   # estimate df that maximizes the chisqlog-likelihood  
   opt = optimize( function(df, y){
      sum(dchisq(y, df=df, log=TRUE))
   }, interval=c(2, 1e6), y=chisq_stats, maximum=TRUE)
   df.approx  = opt$maximum

   # compute p-value based on permuted statistic
   p.value = pchisq( fit$X2, df=df.approx, lower.tail=FALSE)

   # sum(chisq_stats > fit$statistic) / nperm

   list( p.value = as.numeric(p.value),
      statistic = fit$statistic,
      df.approx = df.approx,
      stat_logdet = fit$stat_logdet)
}







