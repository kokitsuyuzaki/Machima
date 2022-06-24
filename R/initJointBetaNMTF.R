# Initialization: W_RNA, H_RNA, W_Epi, H_Epi, T
.initJointBetaNMTF <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        int <- .initJointBetaNMTF_Matrix(X_RNA, X_Epi, T, fixT, pseudocount, J, init)
    }else{
        int <- .initJointBetaNMTF_List(X_RNA, X_Epi, T, fixT, pseudocount, J, init)
    }
    int
}

.initJointBetaNMTF_Matrix <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init){
    X_RNA[which(X_RNA == 0)] <- pseudocount
    X_Epi[which(X_Epi == 0)] <- pseudocount
    if(init == "NMF"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA, J=J), X_RNA)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
    }
    if(init == "Random"){
        W_RNA <- .normalizeCols(
            matrix(runif(nrow(X_RNA)*J),
            nrow=nrow(X_RNA), ncol=J))$A
        H_RNA <- t(W_RNA) %*% X_RNA
    }
    # Random H_Epi/W_Epi
    H_Epi <- .normalizeRows(
        matrix(runif(J*ncol(X_Epi)),
        nrow=J, ncol=ncol(X_Epi)))$A
    W_Epi <- X_Epi %*% t(H_Epi)
    if(!fixT){
        nr <- nrow(X_Epi)
        nc <- ncol(X_RNA)
        T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
    }
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, W_Epi=W_Epi, H_Epi=H_Epi, T=T)
}

.initJointBetaNMTF_List <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init){
    X_RNA <- lapply(X_RNA, function(x){
        x[which(x == 0)] <- pseudocount
        x
    })
    X_Epi <- lapply(X_Epi, function(x){
        x[which(x == 0)] <- pseudocount
        x
    })
    X_RNA2 <- do.call("rbind", X_RNA)
    if(init == "NMF"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA2, J=J), X_RNA2)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA2 <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
    }
    if(init == "Random"){
        W_RNA2 <- .normalizeCols(
            matrix(runif(nrow(X_RNA2)*J),
            nrow=nrow(X_RNA2), ncol=J))$A
        H_RNA <- t(W_RNA2) %*% X_RNA
    }
    # Matrix (W_RNA2) => List (W_RNA)
    dim_RNA <- unlist(lapply(X_RNA, nrow))
    ends_RNA <- cumsum(dim_RNA)
    starts_RNA <- c(1, 1+ends_RNA[1:(length(ends_RNA)-1)])
    W_RNA <- list()
    length(W_RNA) <- length(dim_RNA)
    for(i in seq(dim_RNA)){
        W_RNA[[i]] <- W_RNA2[starts_RNA[i]:ends_RNA[i], ]
    }
    # Random H_Epi/W_Epi
    H_Epi <- .normalizeRows(
            matrix(runif(J*ncol(X_Epi[[1]])),
            nrow=J, ncol=ncol(X_Epi[[1]])))$A
    X_Epi2 <- do.call("rbind", X_Epi)
    W_Epi2 <- X_Epi2 %*% t(H_Epi)
    original_Epi <- unlist(lapply(X_Epi, nrow))
    ends_Epi <- cumsum(original_Epi)
    starts_Epi <- c(1, 1 + ends_Epi[1:(length(ends_Epi)-1)])
    W_Epi <- list()
    length(W_Epi) <- length(original_Epi)
    for(i in seq(original_Epi)){
        W_Epi[[i]] <- W_Epi2[starts_Epi[i]:ends_Epi[i], ]
    }
    if(!fixT){
        T <- lapply(seq_along(X_Epi), function(x){
            nr <- nrow(X_Epi[[x]])
            nc <- nrow(X_RNA[[x]])
            matrix(runif(nr*nc), nrow=nr, ncol=nc)
        })
    }
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, W_Epi=W_Epi, H_Epi=H_Epi, T=T)
}
