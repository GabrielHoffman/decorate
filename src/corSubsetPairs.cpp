// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec corSubsetPairs( const mat & Y, const vec & idxi, const vec & idxj){
  
	// Result is symmatric upper triangle, 
	// so any entries where i > j should be flipped
	// all pairs should be unique
	// but R code is slow
	//df_idx = unique(t(apply(cbind(idxi, idxj), 1, function(x) sort(x))))

	// for (int k=0; k<idxj.n_rows; k++) {
	// 	if( idxi(k) > idxj(k) ){
	// 		std::swap( idxi(k), idxj(k) );
	// 		// Rcpp::Rcout << k << ' '  << idxi(k) << ' ' << idxj(k) << std::endl;
	// 	}
	// }

	// initialize results vector  
	vec rho( idxj.n_rows );

	double N = Y.n_rows;

	for (int k=0; k<idxj.n_rows; k++) {

		// eval correlation for each pair of indecies
		// convert indeces to C from R: so subtract 1
		rho(k) = sum(Y.col( idxi(k) - 1 ) % Y.col( idxj(k) - 1 )) / (N-1);
	}

	return rho; 
}