# Update rules for T in the symmetric model.
# X_Epi ~ G * H_Sym * t(G), G = T * W_RNA.
# dD/dT = dD/dG * t(W_RNA) = 2 * M * G * H_Sym * t(W_RNA)
# The factor 2 cancels since T only appears in the Epi term.

.updateT2 <- function(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root){
    if(is.matrix(X_Epi)){
        T <- .updateT2_Matrix(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root)
    }else{
        T <- .updateT2_List(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root)
    }
    T
}

.updateT2_Matrix <- function(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root){
    G <- T %*% W_RNA
    S_hat <- G %*% H_Sym %*% t(G)
    HW <- H_Sym %*% t(W_RNA)
    numer <- (S_hat^(Beta - 2) * X_Epi) %*% G %*% HW
    if(orthT){
        denom <- T %*% t(T) %*% X_Epi %*% G %*% HW + L1_T + L2_T * T
    }else{
        denom <- S_hat^(Beta - 1) %*% G %*% HW + L1_T + L2_T * T
    }
    T * (numer / denom)^.rho(Beta, root)
}

.updateT2_List <- function(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root){
    lapply(seq_along(X_Epi), function(x){
        G <- T[[x]] %*% W_RNA[[x]]
        S_hat <- G %*% H_Sym %*% t(G)
        HW <- H_Sym %*% t(W_RNA[[x]])
        numer <- (S_hat^(Beta - 2) * X_Epi[[x]]) %*% G %*% HW
        if(orthT){
            denom <- T[[x]] %*% t(T[[x]]) %*% X_Epi[[x]] %*% G %*% HW + L1_T + L2_T * T[[x]]
        }else{
            denom <- S_hat^(Beta - 1) %*% G %*% HW + L1_T + L2_T * T[[x]]
        }
        T[[x]] * (numer / denom)^.rho(Beta, root)
    })
}
