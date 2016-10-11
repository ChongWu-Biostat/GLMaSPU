#' Asymptotic based Sum of Powered Score (SPU) tests and adaptive SPU (aSPU) test.
#'
#' It gives p-values of the SPU tests and aSPU test.
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
#' Chong Wu, Gongjun Xu and Wei Pan, "An Adaptive test on high dimensional parameters in generalized linear models" (Submitted)
#'
#' @examples
#'
#' p = 200
#' n = 100
#' beta = c(1,3,3)
#' s = 0.15
#' non.zero = floor(p * s)
#' signal.r = 0.02
#' seed = 1
#' alpha = c(rep(signal.r,non.zero),rep(0,p-non.zero))
#' dat = generate_data(seed, n = n, p = p, beta = beta,alpha = alpha)
#' cov = dat$Z
#' X = dat$X
#' Y = dat$Y
#' aSPU_apval(Y, X, cov = cov, pow = c(1:6, Inf),resample = "perm", model = "gaussian",  n.perm = 1000)
#' # The p-values are similar to the resample based one
#'
aSPU_apval <- function(Y, X,cov = NULL, pow = c(1:6, Inf), resample = "boot", model= "gaussian",n.perm = 5000){
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    if(is.null(cov)) {
        r <- Y - mean(Y)
        U <- t(X) %*% r
        U <- U / n
        X.new = sweep(X,MARGIN=1,r,`*`)
    } else {
        regression.data <- as.data.frame(cbind(trait = Y, cov))
        if(is.null(colnames(cov))) {
            colnames(regression.data) <- c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(regression.data) <- c("trait", colnames(cov))
        }
        
        fit <- glm(trait~.,family = model, data = regression.data)
        pis <- fitted.values(fit)
        rownames(regression.data) %in% names(pis)
        r <- Y - pis
        X.new <- sweep(X,MARGIN=1,r,`*`)
        U <- t(X) %*% r
        U <- U / n
    }
    
    ## calculate the expecatation, variance, and covariance of L(gamma)
    parametric.boot <- matrix(NA, n.perm, dim(U)[1])
    if(resample == "boot") {
        for (b in 1:n.perm) {
            set.seed(b)
            if (is.null(cov)) {
                Y0 <- sample(Y, length(Y))
                ##  Null score vector:
                parametric.boot[b,]<- t(X) %*% (Y0-mean(Y0))
            } else {
                
                if( model == "gaussian" ) {
                    Y0 <- pis + sample(fit$residuals, n, replace = F )
                    tdat0 <- data.frame(trait=Y0, cov)
                    fit0 <-  glm(trait ~., data = tdat0)
                    yfits0 <- fitted.values(fit0)
                    U0 <- t(X) %*% (Y0 - yfits0)
                    parametric.boot[b,]<- U0/n
                } else {
                    ## with nuisance parameters:
                    Y0 <- rbinom(n,1, prob = pis)
                    # for(i in 1:n) Y0[i] <- sample(c(1,0), 1, prob=c(pis[i], 1-pis[i]) )
                    tdat0<-data.frame(trait=Y0, cov)
                    fit0<-glm(trait~., family=model, data=tdat0)
                    yfits0<-fitted.values(fit0)
                    U0 <- t(X) %*% (Y0 - yfits0)
                    parametric.boot[b,]<- U0/n
                }
            }
        }
    } else {
        for (b in 1:n.perm) {
            r0 <- sample(r, length(r))
            parametric.boot[b,] <-  as.vector(t(X) %*% r0) / n
        }
    }
    
    ##observed statistics
    pval <- numeric(length(pow) + 1)
    L <- numeric(length(pow))
    if(sum(pow==Inf)==0) {
        L.e <- numeric(length(pow))
        L.var <- numeric(length(pow))
        
        stan.L <- numeric(length(pow))
        
        boot.test.stat <- matrix(NA,n.perm,length(pow))
        
        for (b in 1:length(pow)) {
            boot.test.stat[,b] <-rowSums(parametric.boot^{pow[b]})
        }
    } else {
        L.e <- numeric(length(pow) - 1)
        L.var <- numeric(length(pow) - 1)
        
        stan.L <- numeric(length(pow) - 1)
        
        boot.test.stat <- matrix(NA,n.perm,length(pow)-1)
        
        for (b in 1:(length(pow)-1)) {
            boot.test.stat[,b] <-rowSums(parametric.boot^{pow[b]})
        }
    }
    
    L.e <- colMeans(boot.test.stat)
    boot.var <- var(boot.test.stat)
    
    L.var <- diag(boot.var)
    boot.cor <- cor(boot.test.stat)
    
    ########################################
    
    for(i in 1:length(pow)){
        if(pow[i] != Inf){
            L[i] <- sum(U^(pow[i]))
            stan.L[i] <- (L[i] - L.e[i]) / sqrt(L.var[i])
            if(pow[i] %% 2 == 1) pval[i] <- 2 * (1 - pnorm(abs(stan.L[i])))
            if(pow[i] %% 2 == 0) pval[i] <- 1 - pnorm(stan.L[i])
        } else {
            diag.sam.cov <- diag(cov((X.new)))
            diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
            L[i] <- n * max(U^2/diag.sam.cov)
            stan.L.inf <- L[i] - (2*log(p) - log(log(p)))
            pval[i] <- pval.inf <- 1 - exp(-exp(-stan.L.inf/2)/sqrt(pi))
        }
    }
    
    f.pow <- pow[pow != Inf]
    odd.ga <- f.pow[f.pow %% 2 == 1]
    odd.ga.id <- which(f.pow %% 2 == 1)
    even.ga <- f.pow[f.pow %% 2 == 0]
    even.ga.id <- which(f.pow %% 2 == 0)
    n.odd.ga <- length(odd.ga)
    n.even.ga <- length(even.ga)
    R_O <- matrix(NA, n.odd.ga, n.odd.ga)
    R_E <- matrix(NA, n.even.ga, n.even.ga)
    
    R_O <- boot.cor[odd.ga.id,odd.ga.id]
    R_E <- boot.cor[even.ga.id,even.ga.id]
    
    TO <- max(abs(stan.L[odd.ga.id]))
    TE <- max(stan.L[even.ga.id])
    pval_O <- 1 - pmvnorm(lower = -rep(TO, n.odd.ga), upper = rep(TO, n.odd.ga), mean = rep(0, n.odd.ga), sigma = R_O)
    pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(TE, n.even.ga), mean = rep(0, n.even.ga), sigma = R_E)
    
    if(sum(pow==Inf)==0) {
        pval.min <- min(c(pval_O, pval_E))
        pval[length(pow) + 1] <- 1 - (1 - pval.min)^2
        names(pval) <- c(paste("SPU_", pow, sep = ""), "aSPU")
    } else {
        pval.min <- min(c(pval_O, pval_E, pval.inf))
        pval[length(pow) + 1] <- 1 - (1 - pval.min)^3
        names(pval) <- c(paste("SPU_", pow, sep = ""), "aSPU")
    }
    out = list(pvs = pval, Ts = L)
    return(out)
}