.updateW_RNA <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T,
    Beta, L1_W_RNA, L2_W_RNA, orthW_RNA, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        W_RNA <- .updateW_RNA_Matrix(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi,
            T, Beta, L1_W_RNA, L2_W_RNA, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }else{
        W_RNA <- .updateW_RNA_List(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi,
            T, Beta, L1_W_RNA, L2_W_RNA, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }
    W_RNA
}

.updateW_RNA_Matrix <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T,
    Beta, L1_W_RNA, L2_W_RNA, orthW_RNA, root, Pi_RNA, Pi_Epi){
    WH <- W_RNA %*% H_RNA
    numer1 <- Pi_RNA * ((WH^(Beta - 2) * X_RNA) %*% t(H_RNA))
    if(orthW_RNA){
        denom1 <- Pi_RNA * (W_RNA %*% t(W_RNA) %*% X_RNA %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA)
    }else{
        denom1 <- Pi_RNA * ((WH^(Beta - 1) %*% t(H_RNA)) + L2_W_RNA + L2_W_RNA * W_RNA)
    }
    TWH <- T %*% W_RNA %*% H_Epi
    numer2 <- Pi_Epi * (t(T) %*% (TWH^(Beta - 2) * X_Epi) %*% t(H_Epi))
    denom2 <- Pi_Epi * (t(T) %*% TWH^(Beta - 1) %*% t(H_Epi) + L2_W_RNA + L2_W_RNA * W_RNA)
    W_RNA * ((numer1 + numer2)/ (denom1 + denom2))^.rho(Beta, root)
}

.updateW_RNA_List <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T,
    Beta, L1_W_RNA, L2_W_RNA, orthW_RNA, root, Pi_RNA, Pi_Epi){
        lapply(seq_along(X_RNA), function(x){
            WH <- W_RNA[[x]] %*% H_RNA
            numer1 <- Pi_RNA[[x]] * ((WH^(Beta - 2) * X_RNA[[x]]) %*% t(H_RNA))
            if(orthW_RNA){
                denom1 <- Pi_RNA * (W_RNA[[x]] %*% t(W_RNA[[x]]) %*% X_RNA[[x]] %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            }else{
                denom1 <- Pi_RNA * (WH^(Beta - 1) %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            }
            TWH <- T[[x]] %*% W_RNA[[x]] %*% H_Epi
            numer2 <- Pi_Epi[[x]] * (t(T[[x]]) %*% (TWH^(Beta - 2) * X_Epi[[x]]) %*% t(H_Epi))
            denom2 <- Pi_Epi[[x]] * (t(T[[x]]) %*% TWH^(Beta - 1) %*% t(H_Epi) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            W_RNA[[x]] * ((numer1 + numer2)/ (denom1 + denom2))^.rho(Beta, root)
        })
}
