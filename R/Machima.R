#' Title
#'
#' @param X_RNA
#' @param X_Epi
#' @param fixW_RNA
#' @param fixH_RNA
#' @param fixT
#' @param pseudocount
#' @param L1_W_RNA
#' @param L2_W_RNA
#' @param L1_H_RNA
#' @param L2_H_RNA
#' @param L1_T
#' @param L2_T
#' @param L1_H_Epi
#' @param L2_H_Epi
#' @param J
#' @param Beta
#' @param thr
#' @param viz
#' @param figdir
#' @param init
#' @param num.iter
#' @param verbose
#' @return
#' @examples
#'
#' @importFrom nnTensor NMF
#' @importFrom fields image.plot
#' @export
Machima <- function(X_RNA, X_Epi, T=NULL,
    fixW_RNA=TRUE, fixH_RNA=TRUE, fixT=FALSE,
    pseudocount=1e-10,
    L1_W_RNA=1e-10, L2_W_RNA=1e-10,
    L1_H_RNA=1e-10, L2_H_RNA=1e-10,
    L1_T=1e-10, L2_T=1e-10,
    L1_H_Epi=1e-10, L2_H_Epi=1e-10,
    J=3, Beta=2, thr=1e-10, viz=FALSE, figdir=NULL,
    init = c("NMF", "Random"),
    num.iter=30, verbose=FALSE){
    # Argument Check
    init <- match.arg(init)
    .checkJointBetaNMTF(X_RNA, X_Epi, T,
        fixW_RNA, fixH_RNA, fixT, pseudocount,
        L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
        L1_T, L2_T, L1_H_Epi, L2_H_Epi,
        J, Beta, thr, viz, figdir, num.iter, verbose)
    # Initialization
    int <- .initJointBetaNMTF(X_RNA, X_Epi, T, fixT, pseudocount, J, init)
    X_RNA <- int$X_RNA
    X_Epi <- int$X_Epi
    W_RNA <- int$W_RNA
    H_RNA <- int$H_RNA
    W_Epi <- int$W_Epi
    H_Epi <- int$H_Epi
    T <- int$T
    X_GAM <- W_RNA %*% H_Epi
    # Before Update
    if(viz && !is.null(figdir)){
        png(filename = paste0(figdir, "/0.png"),
            width=2500, height=500)
        .multiImagePlots2(list(X_RNA, X_GAM, X_Epi))
        dev.off()
    }
    if(viz && is.null(figdir)){
        .multiImagePlots2(list(X_RNA, X_GAM, X_Epi))
    }
    # Iteration
    iter <- 1
    while (iter <= num.iter){
        # Step1: Update W_RNA
        if(!fixW_RNA){
        W_RNA <- .updateW_RNA(X_RNA, W_RNA, H_RNA, H_Epi, T, Beta, L1_W_RNA, L2_W_RNA)
        }
        # Step2: Update H_RNA
        if(!fixH_RNA){
        H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA)
        }
        # Step3: Update T
        if(!fixT){
            T <- .updateT(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T)
        }
        # Step4: Update H_Epi
        H_Epi <- .updateH_Epi(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi)
        # X_GAM
        X_GAM <- W_RNA %*% H_Epi
        # After Update
        if(viz && !is.null(figdir)){
            png(filename = paste0(figdir, "/", iter, ".png"),
                width=2500, height=500)
            .multiImagePlots2(list(X_RNA, X_GAM, X_Epi))
            dev.off()
        }
        if(viz && is.null(figdir)){
            .multiImagePlots2(list(X_RNA, X_GAM, X_Epi))
        }
        if(verbose){
            cat(paste0(iter, " / ", num.iter, "\n"))
        }
        iter <- iter + 1
    }
    # Output
    list(W_RNA=W_RNA, H_RNA=H_RNA,
        W_Epi=W_Epi, H_Epi=H_Epi,
        X_GAM=X_GAM, T=T)
}

# Check
.checkJointBetaNMTF <- function(X_RNA, X_Epi, T,
    fixW_RNA, fixH_RNA, fixT, pseudocount,
    L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
    L1_T, L2_T, L1_H_Epi, L2_H_Epi,
    J, Beta, thr, viz, figdir, num.iter, verbose){
    # Check X_RNA
    stopifnot(is.matrix(X_RNA))
    # Check X_Epi
    stopifnot(is.matrix(X_Epi))
    # Check T
    check1 <- is.matrix(T)
    check2 <- is.null(T)
    if(!check1 && !check2){
        stop("Please specify T as a matrix or NULL")
    }
    # Check fix
    stopifnot(is.logical(fixW_RNA))
    stopifnot(is.logical(fixH_RNA))
    stopifnot(is.logical(fixT))
    # Check Pseudo-count
    stopifnot(pseudocount >= 0)
    # Check Regularization Parameters
    stopifnot(L1_W_RNA >= 0)
    stopifnot(L2_W_RNA >= 0)
    stopifnot(L1_H_RNA >= 0)
    stopifnot(L2_H_RNA >= 0)
    stopifnot(L1_T >= 0)
    stopifnot(L2_T >= 0)
    stopifnot(L1_H_Epi >= 0)
    stopifnot(L2_H_Epi >= 0)
    # Check J
    stopifnot(J <= min(dim(X_RNA), dim(X_Epi)))
    # Check Beta
    stopifnot(Beta >= 0)
    # Check thr
    stopifnot(thr >= 0)
    # viz
    stopifnot(is.logical(viz))
    # Check figdir
    if(!is.character(figdir) && !is.null(figdir)){
        stop("Please specify the figdir as a string or NULL")
    }
    # Check num.iter
    stopifnot(num.iter >= 0)
    # Check verbose
    stopifnot(is.logical(verbose))
}

# Initialization: W_RNA, H_RNA, W_Epi, H_Epi, T
.initJointBetaNMTF <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init){
    X_RNA[which(X_RNA == 0)] <- pseudocount
    X_Epi[which(X_Epi == 0)] <- pseudocount
    if(init == "NMF"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA, J=J), X_RNA)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
        # NMF with Epigenome
        out1_3 <- .reArrangeOuts(.returnBestNMF(X_Epi, J=J), X_Epi)
        out1_4 <- .normalizeCols(out1_3$V)
        W_Epi <- t(t(out1_3$U) * out1_4$normA)
        H_Epi <- t(out1_4$A)
    }
    if(init == "Random"){
        W_RNA <- .normalizeCols(
            matrix(runif(nrow(X_RNA)*J),
            nrow=nrow(X_RNA), ncol=J))$A
        H_RNA <- t(W_RNA) %*% X_RNA
        H_Epi <- .normalizeRows(
            matrix(runif(J*ncol(X_Epi)),
            nrow=J, ncol=ncol(X_Epi)))$A
        W_Epi <- X_Epi %*% t(H_Epi)
    }
    if(!fixT){
        T <- .returnBestNMFwithV(X=W_Epi,
            initV=t(W_RNA), J=nrow(W_RNA))$U
    }
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, W_Epi=W_Epi, H_Epi=H_Epi, T=T)
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

# Step1: Update W_RNA
.updateW_RNA <- function(X_RNA, W_RNA, H_RNA, H_Epi, T,
    Beta, L1_W_RNA, L2_W_RNA){
    numer1 <- (((W_RNA %*% H_RNA)^(Beta - 2) * X_RNA) %*% t(H_RNA))
    denom1 <- ((W_RNA %*% H_RNA)^(Beta - 1) %*% t(H_RNA)) + L2_W_RNA + L2_W_RNA * W_RNA
    numer2 <- t(T) %*% ((T %*% W_RNA %*% H_Epi)^(Beta - 2) * X_Epi) %*% t(H_Epi)
    denom2 <- t(T) %*% (T %*% W_RNA %*% H_Epi)^(Beta - 1) %*% t(H_Epi) + L2_W_RNA + L2_W_RNA * W_RNA
    W_RNA * ((numer1 / denom1) + (numer2 / denom2))
}


# Step2: Update H_RNA
.updateH_RNA <- function(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA){
    numer1 <- t(t((W_RNA %*% H_RNA)^(Beta - 2) * X_RNA) %*% W_RNA)
    denom1 <- t(t((W_RNA %*% H_RNA)^(Beta - 1)) %*% W_RNA) + L1_H_RNA + L2_H_RNA * H_RNA
    H_RNA * (numer1 / denom1)
}

# Step3: Update T
.updateT <- function(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T){
    numer1 <- ((T %*% W_RNA %*% H_Epi)^(Beta - 2) * X_Epi) %*% t(H_Epi) %*% t(W_RNA)
    denom1 <- (T %*% W_RNA %*% H_Epi)^(Beta - 1) %*% t(H_Epi) %*% t(W_RNA) + L1_T + L2_T * T
    T * numer1 / denom1
}

# Step4: Update H_Epi
.updateH_Epi <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi){
    numer1 <- t(W_RNA) %*% t(T) %*% ((T %*% W_RNA %*% H_Epi)^(Beta - 2) * X_Epi)

    denom1 <- t(W_RNA) %*% t(T) %*% (T %*% W_RNA %*% H_Epi)^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
    H_Epi * numer1 / denom1
}

# Some Helper Functions
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

.multiImagePlots2 <- function(inputList){
    layout(t(1:3))
    image.plot2(inputList[[1]])
    image.plot2(inputList[[2]])
    image.plot2(inputList[[3]])
}

image.plot2 <- function(A, ...){
    image.plot(t(A[nrow(A):1,]), ...)
}
