#' Resample based Sum of Powered Score (SPU) tests and adaptive SPU (aSPU) test.
#'
#' \code{aSPU_perm} returns p-values of the SPU tests and aSPU test.
#'
#' @param Y Response. It can be binary or continuous trait. A vector with length n (number of observations).
#'
#' @param X Genotype or other data; each row for a subject, and each column for a variable of interest. An n by p matrix (n: number of observations, p: number of predictors).
#'
#' @param cov Covariates. An n by q matrix (n: number of observations, q: number of covariates).
#'
#' @param resample Resample methods. "perm" for residual permutations; "boot" for parametric bootstrap.
#'
#' @param model corresponding to the Response. "gaussian" for a quantitative response; "binomial" for a binary response.
#'
#' @param pow Gamma set used in SPU test. A vector of the powers.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @return A list object, Ts : test statistics for the SPU tests and the aSPU test.
#'         pvs : p-values for the SPU and aSPU tests.
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
#' A powerful and adaptive association test for rare variants,
#' Genetics, 197(4), 1081-95
#'
#' @examples
#'
#' p = 200
#' n = 100
#' beta = c(1,3,3)
#' s = 0.15
#' signal.r = 0.08
#' non.zero = floor(p * s)
#' seed = 2
#' alpha = c(rep(signal.r,non.zero),rep(0,p-non.zero))
#' dat = generate_data(seed, n = n, p = p, beta = beta,alpha = alpha)
#' cov = dat$Z
#' X = dat$X
#' Y = dat$Y
#' aSPU_perm(Y, X, cov = cov, pow = c(1:6, Inf),resample = "perm", model = "gaussian",  n.perm = 1000)
#' # The p-values are similar to the asymptotic based one
#'

aSPU_perm <- function(Y, X, cov=NULL, resample = c("perm", "boot"), model=c("gaussian", "binomial"), pow = c(1:6, Inf), n.perm = 1000) {
    X <- as.matrix(X)
    model <- match.arg(model)
    resample <- match.arg(resample)
    
    if(resample == "boot") {
        aSPUboot2_perm(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
    } else {
        aSPUpermC2_perm(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
    }
}
