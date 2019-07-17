
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

# Compute boxM result using only the first *rank* eigen-values 
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

   a = nrow(Y) / ncol(Y) -1
   mats <- aux <- list()
   for(i in seq_len(nlev) ){
      mats[[i]] <-  fxn( Y[group == lev[i], ], method=method )  
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


# boxM_empirical <- function(Y, group, tol=1e-10, nperm = 100){

#    dname <- deparse(substitute(Y))

#    if (!inherits(Y, c("data.frame", "matrix")))
#       stop(paste(dname, "must be a numeric data.frame or matrix!"))
#    if (length(group) != nrow(Y))
#       stop("incompatible dimensions in Y and group!")

#    Y <- as.matrix(Y)
#    gname <- deparse(substitute(group))
#    if (!is.factor(group)) group <- as.factor(as.character(group))

#    valid <- complete.cases(Y, group)

#    if (nrow(Y) > sum(valid)){ 
#       warning(paste(nrow(Y) - sum(valid)), " cases with missing data have been removed.")
#    }
#    Y <- Y[valid,]
#    group <- group[valid]

#    p <- ncol(Y)
#    nlev <- nlevels(group)
#    lev <- levels(group)
#    dfs <- tapply(group, group, length) - 1

#    if( any(dfs < p) ){ 
#       warning("there are one or more levels with less observations than variables!")
#    }

#    # if( useSVD ){

#       # Statistic on real data
#       #------------------------

#       # svd of pooled sample
#       # standardized jointly
#       # Aprime = scale(Y) / sqrt(nrow(Y)-1)

#       compute_m_stat = function(Y, group, tol=1e-10){

#          Y_scale = Y
#          for( key in levels(group)){
#             Y_scale[group==key,] = scale(Y[group==key,])
#          }

#          # standardized by group
#          Aprime = as.matrix(Y_scale / sqrt(nrow(Y_scale)-length(lev)))
#          S_pooled_svd = svd( Aprime, nu=0, nv=0 )   

#          evalues = S_pooled_svd$d^2
#          S_logdet = sum( log(evalues[evalues > tol]) )

#          # svd by group
#          get_Sj_logdet = function(Y, group){
#             Sj_logdet <- vapply( seq_len(nlev), function(i){
#                A = Y[group == lev[i], ]
#                Aprime = scale(A) / sqrt(nrow(A)-1)
#                dcmp = svd( Aprime, nu=0, nv=0)

#                evalues = dcmp$d^2
#                sum( log(evalues[evalues > tol]) )
#             }, numeric(1))      
#             names(Sj_logdet) <- lev
#             Sj_logdet
#          }

#          Sj_logdet <- get_Sj_logdet( Y_scale, group )

#          # compute statistic
#          M <- sum(dfs) * S_logdet - sum(Sj_logdet * dfs)

#          M
#       }

#       M = compute_m_stat(Y, group, tol=tol)

#       # Permutation
#       #-------------

#       # fit statistic for permutated data
#       # get chisq statistic
#       M_perm = vapply( seq_len(nperm), function(i){
#          idx = sample.int(length(group), length(group))
         
#          compute_m_stat(Y, group[idx], tol=tol)
#          }, numeric(1))

#       # estimate df that maximizes the chisq log-likelihood  
#       opt = optimize( function(df, y){
#          sum(dchisq(y, df=df, log=TRUE))
#       }, interval=c(2, 1e6), y=M_perm, maximum=TRUE)
#       df.approx  = opt$maximum

#       # compute p-value based on permuted statistic
#       p.value = pchisq( M, df=df.approx, lower.tail=FALSE)

#       list( p.value   = as.numeric(p.value),
#             statistic = M,
#             df.approx = df.approx)
#    # # }else{

#       # rho = 1
#       # sumabs = 1
#       # trace=FALSE
#       # niter = 50

#       # # sparse PCA
#       #  compute_m_stat = function(Y, group, tol=1e-10){

#       #    Y_scale = Y
#       #    for( key in levels(group)){
#       #       Y_scale[group==key,] = scale(Y[group==key,])
#       #    }

#       #    # standardized by group
#       #    Aprime = as.matrix(Y_scale / sqrt(nrow(Y_scale)-length(lev)))
#       #    S_pooled = crossprod(Aprime)

#       #    dcmp_pooled <- sLED:::symmPMD(S_pooled + rho * diag(p), sumabs = sumabs, trace = trace, niter = niter)

#       #    evalues = dcmp_pooled$d
#       #    S_logdet = sum( log(evalues) )

#       #    # svd by group
#       #    get_Sj_logdet = function(Y, group){
#       #       Sj_logdet <- vapply( seq_len(nlev), function(i){
#       #          A = Y[group == lev[i], ]
#       #          Aprime = scale(A) / sqrt(nrow(A)-1)
#       #          # dcmp = svd( Aprime, nu=0, nv=0)
#       #          Sj = crossprod(Aprime)
#       #          dcmp <- sLED:::symmPMD(Sj + rho * diag(p), sumabs = sumabs, trace = trace, niter = niter)

#       #          evalues = dcmp$d
#       #          sum( log(evalues) )
#       #       }, numeric(1))      
#       #       names(Sj_logdet) <- lev
#       #       Sj_logdet
#       #    }

#       #    Sj_logdet <- get_Sj_logdet( Y_scale, group )

#       #    # compute statistic
#       #    M <- sum(dfs) * S_logdet - sum(Sj_logdet * dfs)

#       #    M
#       # }

#       # M = compute_m_stat(Y, group, tol=tol)

#       # # Permutation
#       # #-------------

#       # # fit statistic for permutated data
#       # # get chisq statistic
#       # M_perm = vapply( seq_len(nperm), function(i){
#       #    idx = sample.int(length(group), length(group))
         
#       #    compute_m_stat(Y, group[idx], tol=tol)
#       #    }, numeric(1))

#       # # estimate df that maximizes the chisq log-likelihood  
#       # opt = optimize( function(df, y){
#       #    sum(dchisq(y, df=df, log=TRUE))
#       # }, interval=c(2, 1e6), y=M_perm, maximum=TRUE)
#       # df.approx  = opt$maximum

#       # # compute p-value based on permuted statistic
#       # p.value = pchisq( M, df=df.approx, lower.tail=FALSE)


#    # }else{
#    #    mats <- aux <- list()
#    #    for(i in seq_len(nlev) ){
#    #       mats[[i]] <- cor(Y[group == lev[i], ])
#    #       aux[[i]] <- mats[[i]] * dfs[i]
#    #    }
#    #    names(mats) <- lev

#    #    S_pooled <- Reduce("+", aux) / sum(dfs)

#    #    # logdet <- log(unlist(lapply(mats, det)))
#    #    Sj_logdet <- vapply(mats, det_partial, tol=tol, numeric(1))

#    #    # logdet_pooled = log(det(pooled))
#    #    S_logdet = det_partial( S_pooled, tol=tol )

#    #    M <- sum(dfs) * S_logdet - sum(Sj_logdet * dfs)

#    # }
# }








