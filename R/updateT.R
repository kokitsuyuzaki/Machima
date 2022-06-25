.updateT <- function(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T, orthT){
    if(is.matrix(X_Epi)){
        T <- .updateT_Matrix(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T, orthT)
    }else{
        T <- .updateT_List(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T, orthT)
    }
    T
}

.updateT_Matrix <- function(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T, orthT){
    TWH <- T %*% W_RNA %*% H_Epi
    HW <- t(H_Epi) %*% t(W_RNA)
    numer1 <- (TWH^(Beta - 2) * X_Epi) %*% HW
    if(orthT){
        denom1 <- T %*% t(T) %*% X_Epi %*% HW + L1_T + L2_T * T
    }else{
        denom1 <- TWH^(Beta - 1) %*% HW + L1_T + L2_T * T
    }
    T * (numer1 / denom1)^.rho(Beta)
}

.updateT_List <- function(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T, orthT){
    lapply(seq_along(X_Epi), function(x){
        TWH <- T[[x]] %*% W_RNA[[x]] %*% H_Epi
        HW <- t(H_Epi) %*% t(W_RNA[[x]])
        numer1 <- (TWH^(Beta - 2) * X_Epi[[x]]) %*% HW
        if(orthT){
            denom1 <- T[[x]] %*% t(T[[x]]) %*% X_Epi[[x]] %*% HW + L1_T + L2_T * T[[x]]
        }else{
            denom1 <- TWH^(Beta - 1) %*% HW + L1_T + L2_T * T[[x]]
        }
        T[[x]] * (numer1 / denom1)^.rho(Beta)
    })
}
