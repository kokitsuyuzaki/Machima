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

.multiImagePlots3 <- function(X_RNA, W_RNA, H_RNA, X_GAM, X_Epi, H_Epi, T){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        recX_RNA <- W_RNA %*% H_RNA
        recX_Epi <- T %*% W_RNA %*% H_Epi
        .multiImagePlots2_Matrix(list(
            X_RNA, recX_RNA, W_RNA, H_RNA, X_GAM,
            X_Epi, recX_Epi, H_Epi, T))
    }else{
        X_RNA2 <- do.call("rbind", X_RNA)
        W_RNA2 <- do.call("rbind", W_RNA)
        if(is.list(X_GAM)){
            X_GAM2 <- do.call("rbind", X_GAM)
        }else{
            X_GAM2 <- X_GAM
        }
        recX_RNA2 <- W_RNA2 %*% H_RNA
        X_Epi2 <- do.call("rbind", X_Epi)
        recX_Epi2 <- do.call("rbind",
            lapply(seq_along(T), function(x){
            T[[x]] %*% W_RNA[[x]] %*% H_Epi
        }))
        .multiImagePlots2_List(list(X_RNA2, recX_RNA2, W_RNA2, H_RNA, X_GAM2,
            X_Epi2, recX_Epi2, H_Epi, T))
    }
}

.multiImagePlots2_Matrix <- function(inputList){
    layout(rbind(1:5, 6:10))
    image.plot2(inputList[[1]], main="X_RNA")
    image.plot2(inputList[[2]], main="rec X_RNA")
    image.plot2(inputList[[3]], main="W_RNA")
    image.plot2(inputList[[4]], main="H_RNA")
    image.plot2(inputList[[5]], main="X_GAM")

    image.plot2(inputList[[6]], main="X_Epi")
    image.plot2(inputList[[7]], main="rec X_Epi")
    image.plot2(inputList[[8]], main="H_Epi")
    image.plot2(inputList[[9]], main="T")
}

.multiImagePlots2_List <- function(inputList){
    l <- length(inputList[[9]])
    m <- length(6:(8+l))
    mylayout <- matrix(0, nrow=2, ncol=m)
    mylayout[1, seq(5)] <- seq(5)
    mylayout[2, ] <- 6:(8+l)
    zeros <- length(which(mylayout == 0))
    layout(mylayout)
    image.plot2(inputList[[1]], main="X_RNA")
    image.plot2(inputList[[2]], main="rec X_RNA")
    image.plot2(inputList[[3]], main="W_RNA")
    image.plot2(inputList[[4]], main="H_RNA")
    image.plot2(inputList[[5]], main="X_GAM")
    .null.plot(zeros)
    image.plot2(inputList[[6]], main="X_Epi")
    image.plot2(inputList[[7]], main="rec X_Epi")
    image.plot2(inputList[[8]], main="H_Epi")
    lapply(seq(l), function(x){
        image.plot2(inputList[[9]][[x]], main=paste0("T_", x))
    })
}

.null.plot <- function(l){
    for(i in seq(l)){
        plot.new()
    }
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

.recErrors <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        d_Beta <- .recErrors_Matrix(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi)
    }else{
        d_Beta <- .recErrors_List(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi)
    }
    d_Beta
}

.recErrors_Matrix <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi){
    left <- Pi_RNA * .BetaDivergence(X_RNA, W_RNA%*%H_RNA, Beta)
    right <- Pi_Epi * .BetaDivergence(X_Epi, T%*%W_RNA%*%H_Epi, Beta)
    left + right
}

.recErrors_List <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi){
    lefts <- sum(unlist(lapply(seq_along(X_RNA), function(x){
        Pi_RNA[[x]] * .BetaDivergence(X_RNA[[x]], W_RNA[[x]]%*%H_RNA, Beta)})))
    rights <- sum(unlist(lapply(seq_along(X_Epi), function(x){
        Pi_Epi[[x]] * .BetaDivergence(X_Epi[[x]], T[[x]]%*%W_RNA[[x]]%*%H_Epi, Beta)})))
    lefts + rights
}

.rho <- function(Beta, root){
    if(root){
        out <- 0.5
    }else{
        if(Beta < 1){
            out <- 1 / (2 - Beta)
        }
        if((1 <= Beta) && (Beta <= 2)){
            out <- 1
        }
        if(Beta > 2){
            out <- 1 / (Beta - 1)
        }
    }
    out
}

.weight <- function(X){
    1 / sum(X^2)
}
