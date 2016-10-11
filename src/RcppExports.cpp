// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// HDGLMC
List HDGLMC(arma::mat X1, arma::mat r1, int nperm);
RcppExport SEXP GLMaSPU_HDGLMC(SEXP X1SEXP, SEXP r1SEXP, SEXP npermSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    __result = Rcpp::wrap(HDGLMC(X1, r1, nperm));
    return __result;
END_RCPP
}
// calcT0C
List calcT0C(arma::mat tXUs, arma::mat yresid, arma::mat powV, int nperm);
RcppExport SEXP GLMaSPU_calcT0C(SEXP tXUsSEXP, SEXP yresidSEXP, SEXP powVSEXP, SEXP npermSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type tXUs(tXUsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type yresid(yresidSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type powV(powVSEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    __result = Rcpp::wrap(calcT0C(tXUs, yresid, powV, nperm));
    return __result;
END_RCPP
}
// GeomanC
List GeomanC(arma::mat X1, arma::mat r1, int nperm);
RcppExport SEXP GLMaSPU_GeomanC(SEXP X1SEXP, SEXP r1SEXP, SEXP npermSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    __result = Rcpp::wrap(GeomanC(X1, r1, nperm));
    return __result;
END_RCPP
}
