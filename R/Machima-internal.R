.Mat2List <- function(X, W2){
    dim <- unlist(lapply(X, nrow))
    ends <- cumsum(dim)
    starts <- c(1, 1+ends[1:(length(ends)-1)])
    W <- list()
    length(W) <- length(dim)
    for(i in seq(dim)){
        W[[i]] <- W2[starts[i]:ends[i], ]
    }
    W
}

.returnBestNMFwithU <- function(X, initU, J){
    outs <- lapply(seq(1), function(x){
        NMF(X, initU=initU, fixU=TRUE, J=J,
            num.iter=30, algorithm="Frobenius")
    })
    bestfit <- unlist(lapply(outs, function(x){rev(x$RecError)[1]}))
    bestfit <- which(bestfit == min(bestfit))[1]
    outs[[bestfit]]
}

.returnBestNMFwithV <- function(X, initV, J){
    outs <- lapply(seq(1), function(x){
        NMF(X, initV=initV, fixV=TRUE, J=J,
            num.iter=30, algorithm="Frobenius")
    })
    bestfit <- unlist(lapply(outs, function(x){rev(x$RecError)[1]}))
    bestfit <- which(bestfit == min(bestfit))[1]
    outs[[bestfit]]
}

.normalizeCols <- function(A){
    normA <- apply(A, 2, function(a){norm(as.matrix(a), "F")})
    A <- t(t(A) / normA)
    list(A=A, normA=normA)
}

.normalizeRows <- function(A){
    normA <- apply(A, 1, function(a){norm(as.matrix(a), "F")})
    A <- A / normA
    list(A=A, normA=normA)
}

.reArrangeOuts <- function(out, X){
    normU <- apply(out$U, 2, function(x){norm(as.matrix(x), "F")})
    normV <- apply(out$V, 2, function(x){norm(as.matrix(x), "F")})
    orderNorm <- order(normU * normV, decreasing=TRUE)
    out$U <- out$U[, orderNorm]
    out$V <- out$V[, orderNorm]
    out
}

.returnBestNMF <- function(X, J){
    outs <- lapply(seq(1), function(x){
        NMF(X, J=J, num.iter=30, algorithm="Frobenius")
    })
    bestfit <- unlist(lapply(outs, function(x){rev(x$RecError)[1]}))
    bestfit <- which(bestfit == min(bestfit))[1]
    outs[[bestfit]]
}

.multiImagePlots3 <- function(X_RNA, X_GAM, X_Epi){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        .multiImagePlots2(list(X_RNA, X_GAM, X_Epi))
    }else{
        X_RNA2 <- do.call("rbind", X_RNA)
        X_Epi2 <- do.call("rbind", X_Epi)
        .multiImagePlots2(list(X_RNA2, X_GAM, X_Epi2))
    }
}

.multiImagePlots2 <- function(inputList){
    layout(t(1:3))
    image.plot2(inputList[[1]])
    image.plot2(inputList[[2]])
    image.plot2(inputList[[3]])
}

image.plot2 <- function(A, ...){
    image.plot(t(A[nrow(A):1,]), ...)
}

.BetaDivergence <- function(X, Y, Beta){
    if(Beta == 1){
        d_Beta <- sum(Y * (log(Y) - log(X)) + (X - Y))
    }else if(Beta == 0){
        d_Beta <- sum((Y / X) - log(Y / X) - 1)
    }else{
        d_Beta <- sum((Y^Beta / Beta * (Beta - 1)) + (X^Beta / Beta) - (Y * X^(Beta-1) / (Beta - 1)))
    }
    d_Beta
}

.recErrors <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        d_Beta <- .recErrors_Matrix(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta)
    }else{
        d_Beta <- .recErrors_List(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta)
    }
    d_Beta
}

.recErrors_Matrix <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta){
    left <- .BetaDivergence(X_RNA, W_RNA%*%H_RNA, Beta)
    right <- .BetaDivergence(X_Epi, T%*%W_RNA%*%H_Epi, Beta)
    left + right
}

.recErrors_List <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta){
    lefts <- sum(unlist(lapply(seq_along(X_RNA), function(x){
        .BetaDivergence(X_RNA[[x]], W_RNA[[x]]%*%H_RNA, Beta)})))
    rights <- sum(unlist(lapply(seq_along(X_Epi), function(x){
        .BetaDivergence(X_Epi[[x]], T[[x]]%*%W_RNA[[x]]%*%H_Epi, Beta)})))
    lefts + rights
}

.rho <- function(Beta){
    if(Beta < 1){
        rho_beta <- 1 / (2 - Beta)
    }
    if((1 <= Beta) && (Beta <= 2)){
        rho_beta <- 1
    }
    if(Beta > 2){
        rho_beta <- 1 / (Beta - 1)
    }
    rho_beta
}
