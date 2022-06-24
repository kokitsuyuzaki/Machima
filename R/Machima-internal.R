
# .returnBestNMFwithU <- function(X, initU, J){
#     outs <- lapply(seq(1), function(x){
#         NMF(X, initU=initU, fixU=TRUE, J=J,
#             num.iter=30, algorithm="Frobenius")
#     })
#     bestfit <- unlist(lapply(outs, function(x){rev(x$RecError)[1]}))
#     bestfit <- which(bestfit == min(bestfit))[1]
#     outs[[bestfit]]
# }

# .returnBestNMFwithV <- function(X, initV, J){
#     outs <- lapply(seq(1), function(x){
#         NMF(X, initV=initV, fixV=TRUE, J=J,
#             num.iter=30, algorithm="Frobenius")
#     })
#     bestfit <- unlist(lapply(outs, function(x){rev(x$RecError)[1]}))
#     bestfit <- which(bestfit == min(bestfit))[1]
#     outs[[bestfit]]
# }

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

.multiImagePlots3 <- function(X_RNA, X_Epi, X_GAM){
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
