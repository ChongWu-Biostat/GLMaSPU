# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// we just rewrite the time consuming part in C.
// [[Rcpp::export()]]
List HDGLMC (arma::mat X1, arma::mat r1,int nperm) {
    //const int n = X1.n_rows;
    const int n = X1.n_cols;
    // containers
    arma::mat T0s(nperm,1);
    T0s.fill(0);
    double T0s1 = 0;
    

    for (int i = 0; i < nperm; i++) {
        arma::mat r10 = shuffle(r1);
        arma::mat U01 = X1 * r10;
        
        T0s1 = accu(pow(U01,2));
        double tmp = 0;
        for (int j = 0; j < n; j++) {
            arma::mat tmp1 = X1.col(j).t() * X1.col(j);
            tmp = tmp + pow(r10(j,0),2) * tmp1(0,0);
        }
        T0s(i,0) = T0s1- tmp;
    }
    List res;
    res["T0s1"] =T0s;
    res["n"] =nperm;
    return(res);
}


// [[Rcpp::export()]]
Rcpp::List calcT0C (arma::mat tXUs, arma::mat yresid, arma::mat diag_cov, arma::mat powV, int nperm) {
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    // containers
    arma::mat T0s1(nperm,npow);
    T0s1.fill(0);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat r10 = shuffle(yresid);
        arma::mat U01 = tXUs * r10;
        
        for (int j = 0; j < npow; j++) {
            if (powV(j,0) == 0) {
                arma::mat tmpU01 = pow(U01,2) / diag_cov;
                T0s1(i,j) = tmpU01.max();
            } else {
                T0s1(i,j) = accu(pow(U01,powV(j,0)));
            }
        }
    }
    Rcpp::List res;
    res["T0"] =T0s1;
    
    return(res);
}


// we just rewrite the time consuming part in C.
// [[Rcpp::export()]]
List GeomanC (arma::mat X1, arma::mat r1, arma::mat diagX1, int nperm) {
    // containers
    arma::mat T0s(nperm,1);
    T0s.fill(0);
    
    
    for (int i = 0; i < nperm; i++) {
        arma::mat r10 = shuffle(r1);
        arma::mat num = r10.t() * X1 * X1.t() * r10;
        
        arma::mat denum = r10.t() * diagX1 * diagX1.t() * r10;
        
        T0s(i,0) = num(0,0) / denum(0,0);
    }
    List res;
    res["T0s1"] =T0s;
    return(res);
}

