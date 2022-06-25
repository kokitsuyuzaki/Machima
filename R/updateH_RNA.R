.updateH_RNA <- function(X_RNA, W_RNA, H_RNA,
            Beta, L1_H_RNA, L2_H_RNA, orthH_RNA){
    if(is.matrix(X_RNA)){
        H_RNA <- .updateH_RNA_Matrix(X_RNA, W_RNA, H_RNA,
            Beta, L1_H_RNA, L2_H_RNA, orthH_RNA)
    }else{
        H_RNA <- .updateH_RNA_List(X_RNA, W_RNA, H_RNA,
            Beta, L1_H_RNA, L2_H_RNA, orthH_RNA)
    }
    H_RNA
}

.updateH_RNA_Matrix <- function(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA, orthH_RNA){
    WH <- W_RNA %*% H_RNA
    numer1 <- t(t(WH^(Beta - 2) * X_RNA) %*% W_RNA)
    if(orthH_RNA){
        denom1 <- t(X_RNA %*% W_RNA %*% H_RNA %*% t(H_RNA)) + L1_H_RNA + L2_H_RNA * H_RNA
    }else{
        denom1 <- t(t(WH^(Beta - 1)) %*% W_RNA) + L1_H_RNA + L2_H_RNA * H_RNA
    }
    H_RNA * (numer1 / denom1)^.rho(Beta)
}

.updateH_RNA_List <- function(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA, orthH_RNA){
    tmp <- lapply(seq_along(X_RNA), function(x){
        WH <- W_RNA[[x]] %*% H_RNA
        numer1 <- t(t(WH^(Beta - 2) * X_RNA[[x]]) %*% W_RNA[[x]])
        if(orthH_RNA){
            denom1 <- t(X_RNA[[x]] %*% W_RNA[[x]] %*% H_RNA %*% t(H_RNA)) + L1_H_RNA + L2_H_RNA * H_RNA
        }else{
            denom1 <- t(t(WH^(Beta - 1)) %*% W_RNA[[x]]) + L1_H_RNA + L2_H_RNA * H_RNA
        }
        (numer1 / denom1)^.rho(Beta)
    })
    H_RNA * Reduce('+', tmp)
}
