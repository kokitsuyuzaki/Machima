# Update rules for H_Sym (J x J, symmetric non-negative).
# Derived from: min D_beta(X_Epi, G * H_Sym * t(G)) where G = T * W_RNA.
# Since G appears on both sides, the gradient w.r.t. H_Sym yields:
#   numer = t(G) %*% (S_hat^(Beta-2) * X_Epi) %*% G
#   denom = t(G) %*% S_hat^(Beta-1) %*% G
# Both numer and denom are symmetric (J x J), preserving H_Sym symmetry.

# --- Tri-factorization mode ---

.updateH_Sym <- function(X_Epi, W_RNA, H_Sym, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_Epi)){
        H_Sym <- .updateH_Sym_Matrix(X_Epi, W_RNA, H_Sym, T, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
    }else{
        H_Sym <- .updateH_Sym_List(X_Epi, W_RNA, H_Sym, T, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
    }
    H_Sym
}

.updateH_Sym_Matrix <- function(X_Epi, W_RNA, H_Sym, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi){
    G <- T %*% W_RNA
    S_hat <- G %*% H_Sym %*% t(G)
    numer <- t(G) %*% (S_hat^(Beta - 2) * X_Epi) %*% G
    if(orthH_Sym){
        A <- t(G) %*% X_Epi %*% G
        denom1 <- (A %*% H_Sym + H_Sym %*% A) / 2 + L1_H_Sym + L2_H_Sym * H_Sym
    }else{
        denom1 <- t(G) %*% S_hat^(Beta - 1) %*% G + L1_H_Sym + L2_H_Sym * H_Sym
    }
    if(orderReg){
        D <- diag(seq(J)) * diag(diag(t(W_RNA) %*% W_RNA))
        denom2 <- (D %*% H_Sym + H_Sym %*% D) / 2
    }else{
        denom2 <- 0
    }
    H_Sym * (numer / (denom1 + denom2))^.rho(Beta, root)
}

.updateH_Sym_List <- function(X_Epi, W_RNA, H_Sym, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi){
    numer <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        G <- T[[x]] %*% W_RNA[[x]]
        S_hat <- G %*% H_Sym %*% t(G)
        Pi_Epi[[x]] * (t(G) %*% (S_hat^(Beta - 2) * X_Epi[[x]]) %*% G)
    }))
    denom1 <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        G <- T[[x]] %*% W_RNA[[x]]
        S_hat <- G %*% H_Sym %*% t(G)
        if(orthH_Sym){
            A <- t(G) %*% X_Epi[[x]] %*% G
            out <- Pi_Epi[[x]] * ((A %*% H_Sym + H_Sym %*% A) / 2 + L1_H_Sym + L2_H_Sym * H_Sym)
        }else{
            out <- Pi_Epi[[x]] * (t(G) %*% S_hat^(Beta - 1) %*% G + L1_H_Sym + L2_H_Sym * H_Sym)
        }
        out
    }))
    denom2 <- Reduce('+',
        lapply(seq_along(X_Epi), function(x){
        if(orderReg){
            D <- diag(seq(J)) * diag(diag(t(W_RNA[[x]]) %*% W_RNA[[x]]))
            out <- Pi_Epi[[x]] * (D %*% H_Sym + H_Sym %*% D) / 2
        }else{
            out <- 0
        }
        out
    }))
    H_Sym * (numer / (denom1 + denom2))^.rho(Beta, root)
}

# --- Horizontal mode ---
# X_GAM = t(T) %*% X_Epi %*% T (n x n, symmetric)
# X_GAM ~ W_RNA %*% H_Sym %*% t(W_RNA)

.updateH_Sym_HZL <- function(X_GAM, W_RNA, H_Sym, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_GAM)){
        H_Sym <- .updateH_Sym_HZL_Matrix(X_GAM, W_RNA, H_Sym, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
    }else{
        H_Sym <- .updateH_Sym_HZL_List(X_GAM, W_RNA, H_Sym, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
    }
    H_Sym
}

.updateH_Sym_HZL_Matrix <- function(X_GAM, W_RNA, H_Sym, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi){
    S_hat <- W_RNA %*% H_Sym %*% t(W_RNA)
    numer <- t(W_RNA) %*% (S_hat^(Beta - 2) * X_GAM) %*% W_RNA
    if(orthH_Sym){
        A <- t(W_RNA) %*% X_GAM %*% W_RNA
        denom1 <- (A %*% H_Sym + H_Sym %*% A) / 2 + L1_H_Sym + L2_H_Sym * H_Sym
    }else{
        denom1 <- t(W_RNA) %*% S_hat^(Beta - 1) %*% W_RNA + L1_H_Sym + L2_H_Sym * H_Sym
    }
    if(orderReg){
        D <- diag(seq(J)) * diag(diag(t(W_RNA) %*% W_RNA))
        denom2 <- (D %*% H_Sym + H_Sym %*% D) / 2
    }else{
        denom2 <- 0
    }
    H_Sym * (numer / (denom1 + denom2))^.rho(Beta, root)
}

.updateH_Sym_HZL_List <- function(X_GAM, W_RNA, H_Sym, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi){
    numer <- Reduce('+',
        lapply(seq_along(X_GAM), function(x){
        S_hat <- W_RNA[[x]] %*% H_Sym %*% t(W_RNA[[x]])
        Pi_Epi[[x]] * (t(W_RNA[[x]]) %*% (S_hat^(Beta - 2) * X_GAM[[x]]) %*% W_RNA[[x]])
    }))
    denom1 <- Reduce('+',
        lapply(seq_along(X_GAM), function(x){
        S_hat <- W_RNA[[x]] %*% H_Sym %*% t(W_RNA[[x]])
        if(orthH_Sym){
            A <- t(W_RNA[[x]]) %*% X_GAM[[x]] %*% W_RNA[[x]]
            out <- Pi_Epi[[x]] * ((A %*% H_Sym + H_Sym %*% A) / 2 + L1_H_Sym + L2_H_Sym * H_Sym)
        }else{
            out <- Pi_Epi[[x]] * (t(W_RNA[[x]]) %*% S_hat^(Beta - 1) %*% W_RNA[[x]] + L1_H_Sym + L2_H_Sym * H_Sym)
        }
        out
    }))
    denom2 <- Reduce('+',
        lapply(seq_along(X_GAM), function(x){
        if(orderReg){
            D <- diag(seq(J)) * diag(diag(t(W_RNA[[x]]) %*% W_RNA[[x]]))
            out <- Pi_Epi[[x]] * (D %*% H_Sym + H_Sym %*% D) / 2
        }else{
            out <- 0
        }
        out
    }))
    H_Sym * (numer / (denom1 + denom2))^.rho(Beta, root)
}
