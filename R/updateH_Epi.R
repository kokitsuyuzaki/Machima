.updateH_Epi <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi){
    if(is.matrix(X_Epi)){
        H_Epi <- .updateH_Epi_Matrix(X_Epi, W_RNA, H_Epi,
            T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi)
    }else{
        H_Epi <- .updateH_Epi_List(X_Epi, W_RNA, H_Epi,
            T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi)
    }
    H_Epi
}

.updateH_Epi_Matrix <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi){
    WT <- t(W_RNA) %*% t(T)
    TWH <- T %*% W_RNA %*% H_Epi
    numer1 <- WT %*% (TWH^(Beta - 2) * X_Epi)
    if(orthH_Epi){
        denom1 <- WT %*% X_Epi %*% t(H_Epi) %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi
    }else{
        denom1 <- WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
    }
    H_Epi * (numer1 / denom1)^.rho(Beta)
}

.updateH_Epi_List <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi){
    tmp <- lapply(seq_along(X_Epi), function(x){
        WT <- t(W_RNA[[x]]) %*% t(T[[x]])
        TWH <- T[[x]] %*% W_RNA[[x]] %*% H_Epi
        numer1 <- WT %*% (TWH^(Beta - 2) * X_Epi[[x]])
        if(orthH_Epi){
            denom1 <- WT %*% X_Epi[[x]] %*% t(H_Epi) %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi
        }else{
            denom1 <- WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
        }
        (numer1 / denom1)^.rho(Beta)
    })
    H_Epi * Reduce('+', tmp)
}
