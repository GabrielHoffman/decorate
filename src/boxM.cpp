// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// rank transform each clumn in Y

arma::mat rankTransform( const arma::mat & Y){

   Function r_rank("rank");   

   NumericVector v;
   arma::vec y;
   arma::uvec idx;
   arma::mat Y_rank(Y.n_rows, Y.n_cols);

   for(int i=0; i<Y.n_cols; i++){
      // get sorted order
      v = r_rank( Y.col(i));

      // save rank to column i
      Y_rank.col(i) = as<arma::vec>(v);
   }

   return Y_rank;
}


// [[Rcpp::export]]
Rcpp::List boxM_fast(const arma::mat & Y, const arma::vec & group, const Rcpp::String & method){

   arma::vec lev = unique(group);
   int nlev = lev.n_elem;

   arma::uvec idx; 
   arma::mat Yi, Ci, C_pooled(Y.n_cols,Y.n_cols);
   arma::vec evalues, dfs(nlev), Si_logDet(nlev), Si_logDet_dfs(nlev); 
   double val, sign;

   for( int i=0; i<nlev; i++){

      idx = find(group == lev(i));
      Yi = Y.rows( idx );

      if( method == "spearman"){
         Ci = cor( rankTransform( Yi ) );
      }else{
         Ci = cor( Yi );
      }

      dfs(i) = idx.n_elem - 1;

      // log det by eigen decomposition
      // evalues = arma::eig_sym( Ci); 
      // Si_logDet(i) = sum( log( evalues ));

      arma::log_det( val, sign, Ci );     
      Si_logDet(i) = val;
      Si_logDet_dfs(i) = Si_logDet(i) * dfs(i);

      if( i == 0 ){
         C_pooled = Ci * dfs(i);
      }else{
         C_pooled = C_pooled + Ci * dfs(i);
      }
   }

   arma::log_det( val, sign, C_pooled / sum(dfs) );   

   // minus2logM <- sum(dfs) * logdet_pooled - sum(logdet * dfs)
   double minus2logM = sum(dfs) * val - sum(Si_logDet_dfs);

   double sum1 = sum(1 / dfs); 
   double p = Y.n_cols;

   double Co = (((2 * p*p) + (3 * p) - 1) / (6 * (p + 1) *
     (nlev - 1))) * (sum1 - (1 / sum(dfs)));

   double X2 = minus2logM * (1 - Co);

   double dfchi = (Rf_choose(p, 2) + p) * (nlev - 1);

   double pval = R::pchisq(X2, dfchi, false, false);

   return Rcpp::List::create( 
                        Named("Si_logDet")  = Si_logDet,
                        Named("dfs")        = dfs,
                        Named("dfchi")      = dfchi,
                        Named("pval")       = pval,
                        Named("X2")         = X2);

}


