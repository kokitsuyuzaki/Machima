.updateH_Epi <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_Epi)){
        H_Epi <- .updateH_Epi_Matrix(X_Epi, W_RNA, H_Epi,
            T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi, root, Pi_RNA, Pi_Epi)
    }else{
        H_Epi <- .updateH_Epi_List(X_Epi, W_RNA, H_Epi,
            T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi, root, Pi_RNA, Pi_Epi)
    }
    H_Epi
}

.updateH_Epi_Matrix <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi, root, Pi_RNA, Pi_Epi){
    WT <- t(W_RNA) %*% t(T)
    TWH <- t(WT) %*% H_Epi
    numer <- WT %*% (TWH^(Beta - 2) * X_Epi)
    if(orthH_Epi){
        denom <- WT %*% X_Epi %*% t(H_Epi) %*% H_Epi + L1_H_Epi + L2_H_Epi * H_Epi
    }else{
        denom <- WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
    }
    H_Epi * (numer / denom)^.rho(Beta, root)
}

.updateH_Epi_List <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi, orthH_Epi, root, Pi_RNA, Pi_Epi){
    numer <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        WT <- t(W_RNA[[x]]) %*% t(T[[x]])
        TWH <- t(WT) %*% H_Epi
        Pi_Epi[[x]] * (WT %*% (TWH^(Beta - 2) * X_Epi[[x]]))
    }))
    denom <- Reduce('+',
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
    H_Epi * (numer / denom)^.rho(Beta, root)
}
