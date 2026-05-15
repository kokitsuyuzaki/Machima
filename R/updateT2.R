# Update rules for T in the symmetric model.
# X_Epi ~ G * H_Sym * t(G), G = T * W_RNA.
# dD/dT = dD/dG * t(W_RNA) = 2 * M * G * H_Sym * t(W_RNA)
# The factor 2 cancels since T only appears in the Epi term.

.updateT2 <- function(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root,
    W_hic=NULL, h_hic=numeric(0), U=NULL){
    if(is.matrix(X_Epi)){
        T <- .updateT2_Matrix(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root,
            W_hic_k=W_hic, h_hic=h_hic, U_k=U)
    }else{
        T <- .updateT2_List(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root,
            W_hic=W_hic, h_hic=h_hic, U=U)
    }
    T
}

.updateT2_Matrix <- function(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root,
    W_hic_k=NULL, h_hic=numeric(0), U_k=NULL){
    W_E <- if(!is.null(U_k)) W_RNA + U_k else W_RNA
    G <- T %*% W_E
    S_hat <- G %*% H_Sym %*% t(G)
    S_hat <- .addHic(S_hat, W_hic_k, h_hic)
    HW <- H_Sym %*% t(W_E)
    numer <- (S_hat^(Beta - 2) * X_Epi) %*% G %*% HW
    if(orthT){
        denom <- T %*% t(T) %*% X_Epi %*% G %*% HW + L1_T + L2_T * T
    }else{
        denom <- S_hat^(Beta - 1) %*% G %*% HW + L1_T + L2_T * T
    }
    T * (numer / denom)^.rho(Beta, root)
}

.updateT2_List <- function(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, L2_T, orthT, root,
    W_hic=NULL, h_hic=numeric(0), U=NULL){
    lapply(seq_along(X_Epi), function(x){
        W_E <- if(!is.null(U)) W_RNA[[x]] + U[[x]] else W_RNA[[x]]
        G <- T[[x]] %*% W_E
        S_hat <- G %*% H_Sym %*% t(G)
        S_hat <- .addHic(S_hat, W_hic[[x]], h_hic)
        HW <- H_Sym %*% t(W_E)
        numer <- (S_hat^(Beta - 2) * X_Epi[[x]]) %*% G %*% HW
        if(orthT){
            denom <- T[[x]] %*% t(T[[x]]) %*% X_Epi[[x]] %*% G %*% HW + L1_T + L2_T * T[[x]]
        }else{
            denom <- S_hat^(Beta - 1) %*% G %*% HW + L1_T + L2_T * T[[x]]
        }
        T[[x]] * (numer / denom)^.rho(Beta, root)
    })
}
