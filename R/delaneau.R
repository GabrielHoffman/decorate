

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


delaneau.test = function( Y, group, method = c("pearson", "kendall", "spearman") ){

      # compute sample level scores
      score = delaneau.score( Y, method )

      # Test of association between score and group variable
      kruskal.test( score, group )   
}


