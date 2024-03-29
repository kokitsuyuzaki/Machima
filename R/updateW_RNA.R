.updateW_RNA <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        W_RNA <- .updateW_RNA_Matrix(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T, J, Beta,
            L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }else{
        W_RNA <- .updateW_RNA_List(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T, J, Beta,
            L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }
    W_RNA
}

.updateW_RNA_Matrix <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
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
    if(orderReg){
        denom3 <- W_RNA %*% (diag(seq(J)) * diag(Pi_RNA * diag(H_RNA %*% t(H_RNA) + Pi_Epi * H_Epi %*% t(H_Epi))))
    }else{
        denom3 <- 0
    }
    W_RNA * ((numer1 + numer2)/ (denom1 + denom2 + denom3))^.rho(Beta, root)
}

.updateW_RNA_List <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
        lapply(seq_along(X_RNA), function(x){
            WH <- W_RNA[[x]] %*% H_RNA
            numer1 <- Pi_RNA[[x]] * ((WH^(Beta - 2) * X_RNA[[x]]) %*% t(H_RNA))
            if(orthW_RNA){
                denom1 <- Pi_RNA[[x]] * (W_RNA[[x]] %*% t(W_RNA[[x]]) %*% X_RNA[[x]] %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            }else{
                denom1 <- Pi_RNA[[x]] * (WH^(Beta - 1) %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            }
            TWH <- T[[x]] %*% W_RNA[[x]] %*% H_Epi
            numer2 <- Pi_Epi[[x]] * (t(T[[x]]) %*% (TWH^(Beta - 2) * X_Epi[[x]]) %*% t(H_Epi))
            denom2 <- Pi_Epi[[x]] * (t(T[[x]]) %*% TWH^(Beta - 1) %*% t(H_Epi) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            if(orderReg){
                denom3 <- W_RNA[[x]] %*% (diag(seq(J)) * diag(Pi_RNA[[x]] * diag(H_RNA %*% t(H_RNA) + Pi_Epi[[x]] * H_Epi %*% t(H_Epi))))
            }else{
                denom3 <- 0
            }
            W_RNA[[x]] * ((numer1 + numer2)/ (denom1 + denom2 + denom3))^.rho(Beta, root)
        })
}

.updateW_RNA_HZL <- function(X_RNA, X_GAM, W_RNA, H_RNA, H_Epi, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        W_RNA <- .updateW_RNA_HZL_Matrix(X_RNA, X_GAM, W_RNA, H_RNA, H_Epi, J, Beta, L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }else{
        W_RNA <- .updateW_RNA_HZL_List(X_RNA, X_GAM, W_RNA, H_RNA, H_Epi, J, Beta, L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }
    W_RNA
}

.updateW_RNA_HZL_Matrix <- function(X_RNA, X_GAM, W_RNA, H_RNA, H_Epi, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    WH1 <- W_RNA %*% H_RNA
    numer1 <- Pi_RNA * ((WH1^(Beta - 2) * X_RNA) %*% t(H_RNA))
    if(orthW_RNA){
        denom1 <- Pi_RNA * (W_RNA %*% t(W_RNA) %*% X_RNA %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA)
    }else{
        denom1 <- Pi_RNA * ((WH1^(Beta - 1) %*% t(H_RNA)) + L2_W_RNA + L2_W_RNA * W_RNA)
    }
    WH2 <- W_RNA %*% H_Epi
    numer2 <- Pi_Epi * ((WH2^(Beta - 2) * X_GAM) %*% t(H_Epi))
    denom2 <- Pi_Epi * (WH2^(Beta - 1) %*% t(H_Epi) + L2_W_RNA + L2_W_RNA * W_RNA)
    if(orderReg){
        denom3 <- W_RNA %*% (diag(seq(J)) * diag(Pi_RNA * diag(H_RNA %*% t(H_RNA) + Pi_Epi * H_Epi %*% t(H_Epi))))
    }else{
        denom3 <- 0
    }
    W_RNA * ((numer1 + numer2)/ (denom1 + denom2 + denom3))^.rho(Beta, root)
}

.updateW_RNA_HZL_List <- function(X_RNA, X_GAM, W_RNA, H_RNA, H_Epi, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
        lapply(seq_along(X_RNA), function(x){
            WH1 <- W_RNA[[x]] %*% H_RNA
            numer1 <- Pi_RNA[[x]] * ((WH1^(Beta - 2) * X_RNA[[x]]) %*% t(H_RNA))
            if(orthW_RNA){
                denom1 <- Pi_RNA[[x]] * (W_RNA[[x]] %*% t(W_RNA[[x]]) %*% X_RNA[[x]] %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            }else{
                denom1 <- Pi_RNA[[x]] * (WH1^(Beta - 1) %*% t(H_RNA) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            }
            WH2 <- W_RNA[[x]] %*% H_Epi
            numer2 <- Pi_Epi[[x]] * ((WH2^(Beta - 2) * X_GAM[[x]]) %*% t(H_Epi))
            denom2 <- Pi_Epi[[x]] * (WH2^(Beta - 1) %*% t(H_Epi) + L2_W_RNA + L2_W_RNA * W_RNA[[x]])
            if(orderReg){
                denom3 <- W_RNA[[x]] %*% (diag(seq(J)) * diag(Pi_RNA[[x]] * diag(H_RNA %*% t(H_RNA) + Pi_Epi[[x]] * H_Epi %*% t(H_Epi))))
            }else{
                denom3 <- 0
            }
            W_RNA[[x]] * ((numer1 + numer2)/ (denom1 + denom2 + denom3))^.rho(Beta, root)
        })
}