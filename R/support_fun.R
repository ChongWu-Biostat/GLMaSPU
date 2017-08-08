aSPUboot2_perm <- function(Y, X, cov=NULL, model=c("gaussian", "binomial"), pow=c(1:8, Inf), n.perm=1000){
    
    model <- match.arg(model)
    
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    if (is.null(cov)){
        ## NO nuisance parameters:
        Xg <- X
        U <- t(Xg) %*% (Y-mean(Y))
        yresids <- Y-mean(Y)
        #        sigma0 = sqrt(sum(yresids^2)/(n-1))
        yfits <- rep(mean(Y), n)
        
    } else {
        ## with nuisance parameters:
        tdat1 <- data.frame(trait=Y, cov)
        
        if(is.null(colnames(cov))) {
            colnames(tdat1) = c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(tdat1) = c("trait", colnames(cov))
        }
        
        fit1 <- glm(trait~., family = model, data=tdat1)
        yfits <- fitted.values(fit1)
        yresids <- Y - yfits
        
        Us <- XUs <- X
        U <- t(XUs) %*% (yresids)
    }
    
    X.new = sweep(X,MARGIN=1,yresids,`*`)
    sam.cov = cov(X.new)
    
    diag.sam.cov <- diag(sam.cov)
    diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
    
    
    
    # test stat's:
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
        Ts[j] = sum(U^pow[j]) else Ts[j] = max(U^2/diag.sam.cov)
    }
    
    # bootstrap:
    T0s = matrix(0, nrow=n.perm, ncol=length(pow))
    Y0 = Y
    for(b in 1:n.perm){
        if (is.null(cov)) {
            Y0 <- sample(Y, length(Y))
            #########Null score vector:
            U0 <- t(Xg) %*% (Y0-mean(Y0))
        } else {
            ## with nuisance parameters:
            if ( model == "gaussian") {
                Y0 <- yfits + sample(yresids, n, replace = F )
                tdat0<-data.frame(trait=Y0, cov)
                fit0<-glm(trait~., data=tdat0)
                yfits0<-fitted.values(fit0)
                U0<-t(XUs) %*% (Y0 - yfits0)
            } else {
                ## with nuisance parameters:
                for(i in 1:n) Y0[i] <- sample(c(1,0), 1, prob=c(yfits[i], 1-yfits[i]) )
                tdat0<-data.frame(trait=Y0, cov)
                fit0<-glm(trait~., family=model, data=tdat0)
                yfits0<-fitted.values(fit0)
                r <- Y0 - yfits0
                U0<-t(XUs) %*% (r)
            }
        }
        
        # test stat's:
        for(j in 1:length(pow))
        if (pow[j] < Inf)
        T0s[b, j] = sum(U0^pow[j]) else T0s[b, j] = max(U0^2/diag.sam.cov)
    }
    
    # bootstrap-based p-values:
    #pPerm0 <- apply( matrix( rep(abs(Ts),n.perm), nrow = n.perm, byrow = T) < abs(T0s), 2, mean)
    pPerm0 = rep(NA,length(pow))
    for ( j in 1:length(pow))
    {
        pPerm0[j] = round( sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm, digits = 8)
        P0s = ( ( n.perm - rank( abs(T0s[,j]) ) ) + 1 ) / (n.perm)
        if (j == 1 ) minp0  = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)
    
    list(Ts = Ts, pvs = pvs)
}


aSPUpermC2_perm <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), pow=c(1:8, Inf), n.perm=1000){
    
    model = match.arg(model)
    
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    if (is.null(cov)){
        ## NO nuisance parameters:
        Xg <- XUs <- X
        U <- t(Xg) %*% (Y-mean(Y))
        yresids <- Y-mean(Y)
        yfits <- rep(mean(Y), n)
        
    } else {
        ## with nuisance parameters:
        tdat1 <- data.frame(trait=Y, cov)
        
        if(is.null(colnames(cov))) {
            colnames(tdat1) = c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(tdat1) = c("trait", colnames(cov))
        }
        
        fit1 <- glm(trait~., family = model, data=tdat1)
        yfits <- fitted.values(fit1)
        yresids <- Y - yfits
        
        
        Us <- XUs <- X
        U <- t(XUs) %*% (yresids)
    }
    X.new = sweep(X,MARGIN=1,yresids,`*`)
    sam.cov = cov(X.new)
    
    diag.sam.cov <- diag(sam.cov)
    diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
    
    
    # test stat's:
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
        Ts[j] = sum(U^pow[j]) else Ts[j] = max(U^2/diag.sam.cov)
    }
    
    pow[pow==Inf] = 0 # pass 0 as infitiy
    T0sC <- calcT0C( t(XUs), as.matrix(yresids),as.matrix(diag.sam.cov), as.matrix(pow), n.perm)
    
    T0s <- T0sC$T0
    pPerm0 = rep(NA,length(pow))
    
    for ( j in 1:length(pow))
    {
        pPerm0[j] = round( sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm, digits = 8)
        P0s = ( ( n.perm - rank( abs(T0s[,j]) ) ) + 1 ) / (n.perm )
        if (j == 1 ) minp0  = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)
    
    list(Ts = Ts, pvs = pvs)
    
}


autocorr.mat <- function(p = 1000, rho = 0.5) {
    mat <- diag(p)
    return(rho^abs(row(mat)-col(mat)))
}

