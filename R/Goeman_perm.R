#' Resample based Goeman test.
#'
#' \code{Goeman_perm} returns resample based p-value for a test proposed by Goeman (2011).
#'
#' \code{Goeman_perm} calculates the resample based p-value. You can calculate the asymptotic based p-value via using R package globaltest. Based on our experience, resample based p-value is often different from the asymptotic based one, except when the dimension of X is larger than the sample size n.
#'
#' @param Y Response. It can be binary or continuous trait. A vector with length n (number of observations).
#'
#' @param X Genotype or other data; each row for a subject, and each column for a variable of interest. An n by p matrix (n: number of observations, p: number of predictors).
#'
#' @param cov Covariates. An n by q matrix (n: number of observations, q: number of covariates).
#'
#' @param model corresponding to the Response. "gaussian" for a quantitative response; "binomial" for a binary response.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @return A list object, Ts : test statistics for the SPU tests and the aSPU test.
#'         pvs : p-values for the SPU and aSPU tests.
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Goeman, J. J., Van Houwelingen, H. C. and Finos, L. (2011). Testing against a high-dimensional
#' alternative in the generalized linear model asymptotic type 1 error control. Biometrika, 98(2), 381-390.
#'
#' @examples
#'
#' p = 200
#' n = 100
#' beta = c(1,3,3)
#' s = 0.15
#' signal.r = 0.02
#' non.zero = floor(p * s)
#' seed = 1
#' alpha = c(rep(signal.r,non.zero),rep(0,p-non.zero))
#' dat = generate_data(seed, n = n, p = p, beta = beta,alpha = alpha)
#' cov = dat$Z
#' X = dat$X
#' Y = dat$Y
#' Goeman_perm(Y, X, cov = cov,model="gaussian", n.perm=1000)
#'
Goeman_perm <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), n.perm=1000){
    
    model = match.arg(model)
    #pow=c(2:8, Inf)
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    #### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        r<-Y-mean(Y)
    } else {
        tdat1 <- data.frame(trait=Y, cov)
        if(is.null(colnames(cov))) {
            colnames(tdat1) = c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(tdat1) = c("trait", colnames(cov))
        }
        fit1 <- glm(trait~.,family=model,data=tdat1)
        pis <- fitted.values(fit1)
        #r<- (Y - pis)/(pis * (1-pis))
        r<- Y - pis
    }
    
    ##observed statistics
    Ts = (t(r) %*% X %*% t(X) %*% r) /(t(r) %*% diag(X) %*% diag(X) %*% r)
    TsC = GeomanC(X,as.matrix(r),n.perm)
    pPerm0=  (sum(TsC$T0s1 >= Ts[1,1])+1) / (n.perm + 1)
    
    res = list(Ts = Ts, pvs = pPerm0)
    res
    
}