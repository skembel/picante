#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat distmat_rcpp(arma::mat comm, arma::mat dis, int threads = 1) {
  
  //Transpose for easy multiplication
  arma::mat tcomm = comm.t();
  
  //Get n communities
  int ncol = tcomm.n_cols;

  //Create distance matrix to allocate
  arma::mat mdistance(ncol, ncol);
  
  // Parallel for speed
#ifdef _OPENMP
  if ( threads > 0 ) {
    omp_set_num_threads( threads );
  }
#endif

  //Loop for distance
#pragma omp parallel for
  for(int i = 0; i < ncol; i++) {
    for(int j = 0; j < ncol; j++) {
      
      //Communities
      arma::vec com_one = tcomm.col(i);
      arma::vec com_two = tcomm.col(j);
      
      //Outer outcome
      arma::mat outer = com_one * com_two.t();
      
      //Multiply and sum accu(1 + 1)
      mdistance(i, j) = accu(outer % dis);
      
    }
  }
  
  return mdistance;
  
}