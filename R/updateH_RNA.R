.updateH_RNA <- function(X_RNA, W_RNA, H_RNA,
            Beta, L1_H_RNA, L2_H_RNA, orthH_RNA, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        H_RNA <- .updateH_RNA_Matrix(X_RNA, W_RNA, H_RNA,
            Beta, L1_H_RNA, L2_H_RNA, orthH_RNA, root, Pi_RNA, Pi_Epi)
    }else{
        H_RNA <- .updateH_RNA_List(X_RNA, W_RNA, H_RNA,
            Beta, L1_H_RNA, L2_H_RNA, orthH_RNA, root, Pi_RNA, Pi_Epi)
    }
    H_RNA
}

.updateH_RNA_Matrix <- function(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA, orthH_RNA, root, Pi_RNA, Pi_Epi){
    WH <- W_RNA %*% H_RNA
    numer <- t(t(WH^(Beta - 2) * X_RNA) %*% W_RNA)
    if(orthH_RNA){
        denom <- (H_RNA %*% t(H_RNA) %*% t(W_RNA) %*% X_RNA) + L1_H_RNA + L2_H_RNA * H_RNA
        t(H_RNA) %*% H_RNA %*% t(X_RNA) %*% W_RNA
    }else{
        denom <- t(t(WH^(Beta - 1)) %*% W_RNA) + L1_H_RNA + L2_H_RNA * H_RNA
    }
    H_RNA * (numer / denom)^.rho(Beta, root)
}

.updateH_RNA_List <- function(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA, orthH_RNA, root, Pi_RNA, Pi_Epi){
    numer <- Reduce('+',
        lapply(seq_along(X_RNA), function(x){
        Pi_RNA[[x]] * (t(t((W_RNA[[x]] %*% H_RNA)^(Beta - 2) * X_RNA[[x]]) %*% W_RNA[[x]]))
    }))
    denom <- Reduce('+',
        lapply(seq_along(X_RNA), function(x){
        if(orthH_RNA){
            out <- Pi_RNA[[x]] * ((H_RNA %*% t(H_RNA) %*% t(W_RNA[[x]]) %*% X_RNA[[x]]) + L1_H_RNA + L2_H_RNA * H_RNA)
        }else{
            out <- Pi_RNA[[x]] * (t(t((W_RNA[[x]] %*% H_RNA)^(Beta - 1)) %*% W_RNA[[x]]) + L1_H_RNA + L2_H_RNA * H_RNA)
        }
        out
    }))
    H_RNA * (numer / denom)^.rho(Beta, root)
}
