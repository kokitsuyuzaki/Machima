.updateH_Epi <- function(X_Epi, W_RNA, H_Epi, T, J, Beta, L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_Epi)){
        H_Epi <- .updateH_Epi_Matrix(X_Epi, W_RNA, H_Epi,
            T, J, Beta, L1_H_Epi, L2_H_Epi, orderReg,
            orthH_Epi, root, Pi_RNA, Pi_Epi)
    }else{
        H_Epi <- .updateH_Epi_List(X_Epi, W_RNA, H_Epi,
            T, J, Beta, L1_H_Epi, L2_H_Epi, orderReg,
            orthH_Epi, root, Pi_RNA, Pi_Epi)
    }
    H_Epi
}

.updateH_Epi_Matrix <- function(X_Epi, W_RNA, H_Epi, T, J, Beta,
    L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi){
    WT <- t(W_RNA) %*% t(T)
    TWH <- t(WT) %*% H_Epi
    numer <- WT %*% (TWH^(Beta - 2) * X_Epi)
    if(orthH_Epi){
        denom1 <- WT %*% X_Epi %*% t(H_Epi) %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi
    }else{
        denom1 <- WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
    }
    if(orderReg){
        denom2 <- (diag(seq(J)) * diag(diag(t(W_RNA) %*% W_RNA))) %*% H_Epi
    }else{
        denom2 <- 0
    }
    H_Epi * (numer / (denom1 + denom2))^.rho(Beta, root)
}

.updateH_Epi_List <- function(X_Epi, W_RNA, H_Epi, T, J, Beta,
    L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi){
    numer <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        WT <- t(W_RNA[[x]]) %*% t(T[[x]])
        TWH <- t(WT) %*% H_Epi
        Pi_Epi[[x]] * (WT %*% (TWH^(Beta - 2) * X_Epi[[x]]))
    }))
    denom1 <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        WT <- t(W_RNA[[x]]) %*% t(T[[x]])
        TWH <- t(WT) %*% H_Epi
        if(orthH_Epi){
            out <- Pi_Epi[[x]] * (WT %*% X_Epi[[x]] %*% t(H_Epi) %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi)
        }else{
            out <- Pi_Epi[[x]] * (WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi)
        }
        out
    }))
    denom2 <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        if(orderReg){
            out <- Pi_Epi[[x]] * (diag(seq(J)) * diag(diag(t(W_RNA[[x]]) %*% W_RNA[[x]]))) %*% H_Epi
        }else{
            out <- 0
        }
        out
    }))
    H_Epi * (numer / (denom1 + denom2))^.rho(Beta, root)
}

.updateH_Epi_HZL <- function(X_GAM, W_RNA, H_Epi, J, Beta, L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_GAM)){
        H_Epi <- .updateH_Epi_HZL_Matrix(X_GAM, W_RNA, H_Epi,
            J, Beta, L1_H_Epi, L2_H_Epi, orderReg,
            orthH_Epi, root, Pi_RNA, Pi_Epi)
    }else{
        H_Epi <- .updateH_Epi_HZL_List(X_GAM, W_RNA, H_Epi,
            J, Beta, L1_H_Epi, L2_H_Epi, orderReg,
            orthH_Epi, root, Pi_RNA, Pi_Epi)
    }
    H_Epi
}

.updateH_Epi_HZL_Matrix <- function(X_GAM, W_RNA, H_Epi, J, Beta, L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi){
    numer <- t(W_RNA) %*% ((W_RNA %*% H_Epi)^(Beta - 2) * X_GAM)
    if(orthH_Epi){
        denom1 <- H_Epi %*% t(X_GAM) %*% W_RNA %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi
    }else{
        denom1 <- t(W_RNA) %*% (W_RNA %*% H_Epi)^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
    }
    if(orderReg){
        denom2 <- (diag(seq(J)) * diag(diag(t(W_RNA) %*% W_RNA))) %*% H_Epi
    }else{
        denom2 <- 0
    }
    H_Epi * (numer / (denom1 + denom2))^.rho(Beta, root)
}

.updateH_Epi_HZL_List <- function(X_GAM, W_RNA, H_Epi, J, Beta, L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi){
    numer <- Reduce('+',
        lapply(seq_along(X_GAM), function(x){
        Pi_Epi[[x]] * (t(W_RNA[[x]]) %*% ((W_RNA[[x]] %*% H_Epi)^(Beta - 2) * X_GAM[[x]]))
    }))
    denom1 <- Reduce('+',
        lapply(seq_along(X_GAM), function(x){
        if(orthH_Epi){
            out <- Pi_Epi[[x]] * (H_Epi %*% t(X_GAM[[x]]) %*% W_RNA[[x]] %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi)
        }else{
            out <- Pi_Epi[[x]] * (t(W_RNA[[x]]) %*% (W_RNA[[x]] %*% H_Epi)^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi)
        }
        out
    }))
    denom2 <- Reduce('+',
        lapply(seq_along(X_GAM), function(x){
        if(orderReg){
            out <- Pi_Epi[[x]] * (diag(seq(J)) * diag(diag(t(W_RNA[[x]]) %*% W_RNA[[x]]))) %*% H_Epi
        }else{
            out <- 0
        }
        out
    }))
    H_Epi * (numer / (denom1 + denom2))^.rho(Beta, root)
}
