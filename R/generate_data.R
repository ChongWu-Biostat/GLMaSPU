#' Generate data for generalized linear models in simulation.
#'
#' \code{generate_data} returns simulated data, including response Y, covariates Z, and variable of interest X.
#'
#' @param seed Random seed.
#'
#' @param n Number of samples
#'
#' @param p Dimension of variable of interest
#'
#' @param beta Coefficients for covariates Z
#'
#' @param alpha Coefficients for variable of interest X
#'
#' @return A list object
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Chong Wu, Gongjun Xu and Wei Pan, "An Adaptive test on high dimensional parameters in generalized linear models" (Submitted)
#'
#' @examples
#'
#' p = 100
#' n = 50
#' beta = c(1,3,3)
#' s = 0.15
#' signal.r = 0.02
#' non.zero = floor(p * s)
#' seed = 1
#' alpha = c(rep(signal.r,non.zero),rep(0,p-non.zero))
#' dat = generate_data(seed, n = n, p = p, beta = beta,alpha = alpha)
#' #X, Y, cov
#' #dat$X; dat$Y; dat$cov
#'
generate_data <- function(seed,n,p,beta, alpha) {
    set.seed(seed)
    Z = cbind(rnorm(n,  0, 1),rnorm(n,  0, 1))
    Z = Z- rep(1, nrow(Z)) %*% t(colMeans(Z))

    sigma = diag(x = sqrt(2),p,p) %*% autocorr.mat(p = p,rho = 0.4) %*% diag(x = sqrt(2),p,p)
    X = mvrnorm(n= n,mu = rep(0,p),Sigma = sigma)
    X = X - rep(1, nrow(X)) %*% t(colMeans(X))
    
    intercept = beta[1]
    beta1 = beta[-1]
    error = rnorm(n, 0,0.5)
    Y = intercept + Z%*% beta1 + X%*% alpha + error
    
    true.cov = 0.01 * sigma
    out = list(Z= Z, X = X, Y = Y, true.cov = true.cov)
    out
}
