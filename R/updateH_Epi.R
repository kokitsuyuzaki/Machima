# Step4: Update H_Epi
.updateH_Epi <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        H_Epi <- .updateH_Epi_Matrix(X_Epi, W_RNA, H_Epi,
            T, Beta, L1_H_Epi, L2_H_Epi)
    }else{
        H_Epi <- .updateH_Epi_List(X_Epi, W_RNA, H_Epi,
            T, Beta, L1_H_Epi, L2_H_Epi)
    }
    H_Epi
}

.updateH_Epi_Matrix <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi){
    WT <- t(W_RNA) %*% t(T)
    TWH <- T %*% W_RNA %*% H_Epi
    numer1 <- WT %*% (TWH^(Beta - 2) * X_Epi)
    denom1 <- WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
    H_Epi * numer1 / denom1
}

.updateH_Epi_List <- function(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi){
    tmp <- lapply(seq_along(X_Epi), function(x){
        WT <- t(W_RNA[[x]]) %*% t(T[[x]])
        TWH <- T[[x]] %*% W_RNA[[x]] %*% H_Epi
        numer1 <- WT %*% (TWH^(Beta - 2) * X_Epi[[x]])
        denom1 <- WT %*% TWH^(Beta - 1) + L1_H_Epi + L2_H_Epi * H_Epi
        numer1 / denom1
    })
    H_Epi * Reduce('+', tmp)
}
