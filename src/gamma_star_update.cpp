#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec gamma_star_update(arma::vec gamma,
                            arma::vec delta1_old,
                            double A21_old,
                            double A22_old,
                            arma::vec delta2_old){

arma::vec gamma_star(gamma.size()); gamma_star.fill(0.00);
arma::vec eta = A21_old*delta1_old +
                A22_old*delta2_old;
int p_z = gamma.size();

for(int j = 0; j < p_z; ++j){
   
   if(gamma(j) == 1.00){
     gamma_star(j) = rnorm_trunc(eta(j),
                                 1.00,
                                 0.00,
                                 datum::inf);
     }
   
   if(gamma(j) == 0.00){
     gamma_star(j) = rnorm_trunc(eta(j),
                                 1.00,
                                 -datum::inf,
                                 0.00);
     }
   
   }
   
return(gamma_star);

}



