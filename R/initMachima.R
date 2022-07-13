.initMachima <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        int <- .initMachima_Matrix(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr)
    }else{
        int <- .initMachima_List(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr)
    }
    int
}

.initMachima_Matrix <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr){
    X_RNA[which(X_RNA == 0)] <- pseudocount
    X_Epi[which(X_Epi == 0)] <- pseudocount
    if(init == "RandomEpi"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA, J=J), X_RNA)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
        # Random H_Epi
        H_Epi <- .normalizeRows(
            matrix(runif(J*ncol(X_Epi)),
            nrow=J, ncol=ncol(X_Epi)))$A
        if(!fixT){
            nr <- nrow(X_Epi)
            nc <- nrow(X_RNA)
            T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
            # # NMF with X_Epi
            # WH <- W_RNA %*% H_Epi
            # T <- .returnBestNMFwithV(X=X_Epi, initV=t(WH), J=nrow(X_RNA))$U
        }
    }
    if(init == "RandomRNA"){
        # Random W_RNA/H_RNA
        W_RNA <- .normalizeCols(
            matrix(runif(nrow(X_RNA)*J),
            nrow=nrow(X_RNA), ncol=J))$A
        H_RNA <- t(W_RNA) %*% X_RNA
        # Random H_Epi
        H_Epi <- .normalizeRows(
            matrix(runif(J*ncol(X_Epi)),
            nrow=J, ncol=ncol(X_Epi)))$A
        # Random T
        if(!fixT){
            nr <- nrow(X_Epi)
            nc <- nrow(X_RNA)
            T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
        }
    }
    if(init == "Random"){
        # Random W_RNA/H_RNA
        W_RNA <- .normalizeCols(
            matrix(runif(nrow(X_RNA)*J),
            nrow=nrow(X_RNA), ncol=J))$A
        H_RNA <- t(W_RNA) %*% X_RNA
        # Random H_Epi
        H_Epi <- .normalizeRows(
            matrix(runif(J*ncol(X_Epi)),
            nrow=J, ncol=ncol(X_Epi)))$A
        # Random T
        if(!fixT){
            nr <- nrow(X_Epi)
            nc <- nrow(X_RNA)
            T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
        }
    }
    if(init == "NMFAlign"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA, J=J), X_RNA)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
        # NMF with Epigenome
        out1_3 <- .reArrangeOuts(.returnBestNMF(X_Epi, J=J), X_Epi)
        out1_4 <- .normalizeCols(out1_3$V)
        TW_Epi <- t(t(out1_3$U) * out1_4$normA)
        H_Epi <- t(out1_4$A)
        # NMF with TW_Epi
        if(!fixT){
            T <- .returnBestNMFwithV(X=TW_Epi, initV=t(W_RNA), J=nrow(X_RNA))$U
        }
    }
    if(init == "NMFAlign2"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA, J=J), X_RNA)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
        # NMF with Epigenome
        out2_1 <- .returnBestNMF(X_Epi, J=J)
        out2_2 <- .normalizeCols(out2_1$V)
        TW_Epi <- t(t(out2_1$U) * out2_2$normA)
        H_Epi <- t(out2_2$A)
        # NMF with TW_Epi
        out2_3 <- .returnBestNMF(TW_Epi, J=nrow(X_RNA))
        out2_4 <- .normalizeCols(out2_3$U)
        W_Epi <- t(out2_3$V * out2_4$normA)
        # Matching by Correlation Coefficient
        cor.matrix <- cor(W_RNA, W_Epi)
        cor.matrix[which(is.na(cor.matrix))] <- 0
        rownames(cor.matrix) <- paste0("W_RNA", seq(J))
        colnames(cor.matrix) <- paste0("W_Epi", seq(J))
        # Alignment
        g <- graph_from_incidence_matrix(cor.matrix, weighted=TRUE)
        W_Epi_idx <- as.vector(max_bipartite_match(g)$matching[paste0("W_RNA", seq(J))])
        W_Epi_idx[which(is.na(W_Epi_idx))] <- setdiff(
            paste0("W_Epi", seq(J)),
            W_Epi_idx[which(!is.na(W_Epi_idx))])
        idx <- unlist(lapply(W_Epi_idx, function(x){
            which(colnames(cor.matrix) == x)
        }))
        W_Epi <- W_Epi[, idx]
        H_Epi <- H_Epi[idx, ]
        WH_Epi <- t(W_Epi %*% H_Epi)
        # NMF with W_Epi
        if(!fixT){
            T <- .returnBestNMFwithV(X=X_Epi, initV=WH_Epi, J=nrow(X_RNA))$U
        }
    }
    # Weight
    Pi_RNA <- .weight(X_RNA)
    Pi_Epi <- .weight(X_Epi)
    # Error
    RecError <- c()
    RelChange <- c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, H_Epi=H_Epi,
        T=T, Pi_RNA=Pi_RNA, Pi_Epi=Pi_Epi,
        RecError=RecError, RelChange=RelChange)
}

.initMachima_List <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr){
    X_RNA <- lapply(X_RNA, function(x){
        x[which(x == 0)] <- pseudocount
        x
    })
    X_Epi <- lapply(X_Epi, function(x){
        x[which(x == 0)] <- pseudocount
        x
    })
    X_RNA2 <- do.call("rbind", X_RNA)
    X_Epi2 <- do.call("rbind", X_Epi)
    if(init == "RandomEpi"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA2, J=J), X_RNA2)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- .Mat2List(X_RNA, out1_1$U)
        H_RNA <- t(out1_2$A * out1_2$normA)
        # Random H_Epi
        H_Epi <- .normalizeRows(
                matrix(runif(J*ncol(X_Epi[[1]])),
                nrow=J, ncol=ncol(X_Epi[[1]])))$A
        # Random T
        if(!fixT){
            T <- lapply(seq_along(X_Epi), function(x){
                nr <- nrow(X_Epi[[x]])
                nc <- nrow(X_RNA[[x]])
                matrix(runif(nr*nc), nrow=nr, ncol=nc)
            })
        }
    }
    if(init == "RandomRNA"){
        # Random W_RNA/H_RNA
        W_RNA2 <- matrix(runif(nrow(X_RNA2)*J), nrow=nrow(X_RNA2), ncol=J)
        W_RNA <- .Mat2List(X_RNA, W_RNA2)
        H_RNA <- t(W_RNA2) %*% X_RNA2
        # NMF with Epigenome
        out1_3 <- .reArrangeOuts(.returnBestNMF(X_Epi2, J=J), X_Epi2)
        out1_4 <- .normalizeCols(out1_3$V)
        W_Epi <- .Mat2List(X_Epi, t(t(out1_3$U) * out1_4$normA))
        H_Epi <- t(out1_4$A)
        # Random T
        if(!fixT){
            T <- lapply(seq_along(X_Epi), function(x){
                nr <- nrow(X_Epi[[x]])
                nc <- nrow(X_RNA[[x]])
                matrix(runif(nr*nc), nrow=nr, ncol=nc)
            })
        }
    }
    if(init == "Random"){
        # Random W_RNA/H_RNA
        W_RNA2 <- matrix(runif(nrow(X_RNA2)*J), nrow=nrow(X_RNA2), ncol=J)
        W_RNA <- .Mat2List(X_RNA, W_RNA2)
        H_RNA <- t(W_RNA2) %*% X_RNA2
        # Random H_Epi
        H_Epi <- .normalizeRows(
                matrix(runif(J*ncol(X_Epi[[1]])),
                nrow=J, ncol=ncol(X_Epi[[1]])))$A
        # Random T
        if(!fixT){
            T <- lapply(seq_along(X_Epi), function(x){
                nr <- nrow(X_Epi[[x]])
                nc <- nrow(X_RNA[[x]])
                matrix(runif(nr*nc), nrow=nr, ncol=nc)
            })
        }
    }
    if(init == "NMFAlign"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA2, J=J), X_RNA2)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- .Mat2List(X_RNA, out1_1$U)
        H_RNA <- t(out1_2$A * out1_2$normA)
        # NMF with Epigenome
        out1_3 <- .reArrangeOuts(.returnBestNMF(X_Epi2, J=J), X_Epi2)
        out1_4 <- .normalizeCols(out1_3$V)
        W_Epi <- .Mat2List(X_Epi, t(t(out1_3$U) * out1_4$normA))
        H_Epi <- t(out1_4$A)
        # NMF with W_Epi
        if(!fixT){
            T <- lapply(seq_along(W_Epi), function(x){
                .returnBestNMFwithV(X=W_Epi[[x]], initV=t(W_RNA[[x]]), J=nrow(X_RNA[[x]]))$U
            })
        }
    }
    # Weight
    Pi_RNA <- lapply(X_RNA, .weight)
    Pi_Epi <- lapply(X_Epi, .weight)
    # Error
    RecError <- c()
    RelChange <- c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, H_Epi=H_Epi,
        T=T, Pi_RNA=Pi_RNA, Pi_Epi=Pi_Epi,
        RecError=RecError, RelChange=RelChange)
}
