# Constrained MU for diagonal H_Sym (sym-CP / sym-PARAFAC).
# Updates h_i directly: h_i <- h_i * (num_i / denom_i)^rho
# where num_i = g_i^T (S^(b-2) * X) g_i, denom_i = g_i^T S^(b-1) g_i

# --- Tri-factorization mode ---

.updateH_Sym_diag <- function(X_Epi, W_RNA, h, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
    W_hic=NULL, h_hic=numeric(0)){
    if(is.matrix(X_Epi)){
        h <- .updateH_Sym_diag_Matrix(X_Epi, W_RNA, h, T, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
            W_hic_k=W_hic, h_hic=h_hic)
    }else{
        h <- .updateH_Sym_diag_List(X_Epi, W_RNA, h, T, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
            W_hic=W_hic, h_hic=h_hic)
    }
    h
}

.updateH_Sym_diag_Matrix <- function(X_Epi, W_RNA, h, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
    W_hic_k=NULL, h_hic=numeric(0)){
    G <- T %*% W_RNA
    S_hat <- G %*% diag(h, nrow=J) %*% t(G)
    S_hat <- .addHic(S_hat, W_hic_k, h_hic)
    numer <- rep(0, J)
    denom <- rep(0, J)
    SbX <- S_hat^(Beta - 2) * X_Epi
    Sb <- S_hat^(Beta - 1)
    for(i in seq_len(J)){
        gi <- G[, i]
        numer[i] <- as.numeric(t(gi) %*% SbX %*% gi)
        denom[i] <- as.numeric(t(gi) %*% Sb %*% gi) + L1_H_Sym + L2_H_Sym * h[i]
    }
    if(orderReg){
        WtW_diag <- diag(t(W_RNA) %*% W_RNA)
        denom <- denom + seq(J) * WtW_diag * h
    }
    h * (numer / denom)^.rho(Beta, root)
}

.updateH_Sym_diag_List <- function(X_Epi, W_RNA, h, T, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
    W_hic=NULL, h_hic=numeric(0)){
    numer <- rep(0, J)
    denom <- rep(0, J)
    for(k in seq_along(X_Epi)){
        G <- T[[k]] %*% W_RNA[[k]]
        S_hat <- G %*% diag(h, nrow=J) %*% t(G)
        S_hat <- .addHic(S_hat, W_hic[[k]], h_hic)
        SbX <- S_hat^(Beta - 2) * X_Epi[[k]]
        Sb <- S_hat^(Beta - 1)
        for(i in seq_len(J)){
            gi <- G[, i]
            numer[i] <- numer[i] + Pi_Epi[[k]] * as.numeric(t(gi) %*% SbX %*% gi)
            denom[i] <- denom[i] + Pi_Epi[[k]] * (as.numeric(t(gi) %*% Sb %*% gi) + L1_H_Sym + L2_H_Sym * h[i])
        }
    }
    if(orderReg){
        for(k in seq_along(W_RNA)){
            WtW_diag <- diag(t(W_RNA[[k]]) %*% W_RNA[[k]])
            denom <- denom + Pi_Epi[[k]] * seq(J) * WtW_diag * h
        }
    }
    h * (numer / denom)^.rho(Beta, root)
}

# --- Horizontal mode ---

.updateH_Sym_diag_HZL <- function(X_GAM, W_RNA, h, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
    W_hic=NULL, h_hic=numeric(0)){
    if(is.matrix(X_GAM)){
        h <- .updateH_Sym_diag_HZL_Matrix(X_GAM, W_RNA, h, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
            W_hic_k=W_hic, h_hic=h_hic)
    }else{
        h <- .updateH_Sym_diag_HZL_List(X_GAM, W_RNA, h, J, Beta,
            L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
            W_hic=W_hic, h_hic=h_hic)
    }
    h
}

.updateH_Sym_diag_HZL_Matrix <- function(X_GAM, W_RNA, h, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
    W_hic_k=NULL, h_hic=numeric(0)){
    S_hat <- W_RNA %*% diag(h, nrow=J) %*% t(W_RNA)
    S_hat <- .addHic(S_hat, W_hic_k, h_hic)
    numer <- rep(0, J)
    denom <- rep(0, J)
    SbX <- S_hat^(Beta - 2) * X_GAM
    Sb <- S_hat^(Beta - 1)
    for(i in seq_len(J)){
        wi <- W_RNA[, i]
        numer[i] <- as.numeric(t(wi) %*% SbX %*% wi)
        denom[i] <- as.numeric(t(wi) %*% Sb %*% wi) + L1_H_Sym + L2_H_Sym * h[i]
    }
    if(orderReg){
        WtW_diag <- diag(t(W_RNA) %*% W_RNA)
        denom <- denom + seq(J) * WtW_diag * h
    }
    h * (numer / denom)^.rho(Beta, root)
}

.updateH_Sym_diag_HZL_List <- function(X_GAM, W_RNA, h, J, Beta,
    L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
    W_hic=NULL, h_hic=numeric(0)){
    numer <- rep(0, J)
    denom <- rep(0, J)
    for(k in seq_along(X_GAM)){
        S_hat <- W_RNA[[k]] %*% diag(h, nrow=J) %*% t(W_RNA[[k]])
        S_hat <- .addHic(S_hat, W_hic[[k]], h_hic)
        SbX <- S_hat^(Beta - 2) * X_GAM[[k]]
        Sb <- S_hat^(Beta - 1)
        for(i in seq_len(J)){
            wi <- W_RNA[[k]][, i]
            numer[i] <- numer[i] + Pi_Epi[[k]] * as.numeric(t(wi) %*% SbX %*% wi)
            denom[i] <- denom[i] + Pi_Epi[[k]] * (as.numeric(t(wi) %*% Sb %*% wi) + L1_H_Sym + L2_H_Sym * h[i])
        }
    }
    if(orderReg){
        for(k in seq_along(W_RNA)){
            WtW_diag <- diag(t(W_RNA[[k]]) %*% W_RNA[[k]])
            denom <- denom + Pi_Epi[[k]] * seq(J) * WtW_diag * h
        }
    }
    h * (numer / denom)^.rho(Beta, root)
}
