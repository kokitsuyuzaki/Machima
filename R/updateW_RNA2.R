# Update rules for W_RNA in the symmetric model.
# X_Epi ~ G * H_Sym * t(G), G = T * W_RNA.
# Since W_RNA appears on both sides through G, the gradient has a factor of 2:
#   dD/dW = dD/dG * dG/dW = 2 * t(T) * M * G * H_Sym
# where M = S_hat^(Beta-1) - (S_hat^(Beta-2) * X_Epi).
# The factor 2 is significant because it weights the Epi contribution
# relative to the RNA term in the combined update.

# --- Tri-factorization mode ---

.updateW_RNA2 <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        W_RNA <- .updateW_RNA2_Matrix(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
            L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }else{
        W_RNA <- .updateW_RNA2_List(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
            L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }
    W_RNA
}

.updateW_RNA2_Matrix <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    # RNA term (standard NMF)
    WH <- W_RNA %*% H_RNA
    numer1 <- Pi_RNA * ((WH^(Beta - 2) * X_RNA) %*% t(H_RNA))
    if(orthW_RNA){
        denom1 <- Pi_RNA * (W_RNA %*% t(W_RNA) %*% X_RNA %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA)
    }else{
        denom1 <- Pi_RNA * ((WH^(Beta - 1) %*% t(H_RNA)) + L1_W_RNA + L2_W_RNA * W_RNA)
    }
    # Epi term (symmetric, factor 2)
    G <- T %*% W_RNA
    GH <- G %*% H_Sym
    S_hat <- GH %*% t(G)
    numer2 <- 2 * Pi_Epi * (t(T) %*% (S_hat^(Beta - 2) * X_Epi) %*% GH)
    denom2 <- 2 * Pi_Epi * (t(T) %*% S_hat^(Beta - 1) %*% GH + L1_W_RNA + L2_W_RNA * W_RNA)
    if(orderReg){
        denom3 <- W_RNA %*% (diag(seq(J)) * diag(Pi_RNA * diag(H_RNA %*% t(H_RNA))))
    }else{
        denom3 <- 0
    }
    W_RNA * ((numer1 + numer2) / (denom1 + denom2 + denom3))^.rho(Beta, root)
}

.updateW_RNA2_List <- function(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    lapply(seq_along(X_RNA), function(x){
        # RNA term
        WH <- W_RNA[[x]] %*% H_RNA
        numer1 <- Pi_RNA[[x]] * ((WH^(Beta - 2) * X_RNA[[x]]) %*% t(H_RNA))
        if(orthW_RNA){
            denom1 <- Pi_RNA[[x]] * (W_RNA[[x]] %*% t(W_RNA[[x]]) %*% X_RNA[[x]] %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA[[x]])
        }else{
            denom1 <- Pi_RNA[[x]] * (WH^(Beta - 1) %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA[[x]])
        }
        # Epi term (symmetric, factor 2)
        G <- T[[x]] %*% W_RNA[[x]]
        GH <- G %*% H_Sym
        S_hat <- GH %*% t(G)
        numer2 <- 2 * Pi_Epi[[x]] * (t(T[[x]]) %*% (S_hat^(Beta - 2) * X_Epi[[x]]) %*% GH)
        denom2 <- 2 * Pi_Epi[[x]] * (t(T[[x]]) %*% S_hat^(Beta - 1) %*% GH + L1_W_RNA + L2_W_RNA * W_RNA[[x]])
        if(orderReg){
            denom3 <- W_RNA[[x]] %*% (diag(seq(J)) * diag(Pi_RNA[[x]] * diag(H_RNA %*% t(H_RNA))))
        }else{
            denom3 <- 0
        }
        W_RNA[[x]] * ((numer1 + numer2) / (denom1 + denom2 + denom3))^.rho(Beta, root)
    })
}

# --- Horizontal mode ---
# X_GAM = t(T) %*% X_Epi %*% T (n x n, symmetric)
# X_GAM ~ W_RNA %*% H_Sym %*% t(W_RNA)
# W_RNA appears on both sides, so factor 2 applies.

.updateW_RNA2_HZL <- function(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        W_RNA <- .updateW_RNA2_HZL_Matrix(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
            L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }else{
        W_RNA <- .updateW_RNA2_HZL_List(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
            L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
    }
    W_RNA
}

.updateW_RNA2_HZL_Matrix <- function(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    # RNA term
    WH1 <- W_RNA %*% H_RNA
    numer1 <- Pi_RNA * ((WH1^(Beta - 2) * X_RNA) %*% t(H_RNA))
    if(orthW_RNA){
        denom1 <- Pi_RNA * (W_RNA %*% t(W_RNA) %*% X_RNA %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA)
    }else{
        denom1 <- Pi_RNA * ((WH1^(Beta - 1) %*% t(H_RNA)) + L1_W_RNA + L2_W_RNA * W_RNA)
    }
    # Epi term (symmetric, factor 2)
    WH2 <- W_RNA %*% H_Sym
    S_hat <- WH2 %*% t(W_RNA)
    numer2 <- 2 * Pi_Epi * ((S_hat^(Beta - 2) * X_GAM) %*% WH2)
    denom2 <- 2 * Pi_Epi * (S_hat^(Beta - 1) %*% WH2 + L1_W_RNA + L2_W_RNA * W_RNA)
    if(orderReg){
        denom3 <- W_RNA %*% (diag(seq(J)) * diag(Pi_RNA * diag(H_RNA %*% t(H_RNA))))
    }else{
        denom3 <- 0
    }
    W_RNA * ((numer1 + numer2) / (denom1 + denom2 + denom3))^.rho(Beta, root)
}

.updateW_RNA2_HZL_List <- function(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi){
    lapply(seq_along(X_RNA), function(x){
        # RNA term
        WH1 <- W_RNA[[x]] %*% H_RNA
        numer1 <- Pi_RNA[[x]] * ((WH1^(Beta - 2) * X_RNA[[x]]) %*% t(H_RNA))
        if(orthW_RNA){
            denom1 <- Pi_RNA[[x]] * (W_RNA[[x]] %*% t(W_RNA[[x]]) %*% X_RNA[[x]] %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA[[x]])
        }else{
            denom1 <- Pi_RNA[[x]] * (WH1^(Beta - 1) %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA[[x]])
        }
        # Epi term (symmetric, factor 2)
        WH2 <- W_RNA[[x]] %*% H_Sym
        S_hat <- WH2 %*% t(W_RNA[[x]])
        numer2 <- 2 * Pi_Epi[[x]] * ((S_hat^(Beta - 2) * X_GAM[[x]]) %*% WH2)
        denom2 <- 2 * Pi_Epi[[x]] * (S_hat^(Beta - 1) %*% WH2 + L1_W_RNA + L2_W_RNA * W_RNA[[x]])
        if(orderReg){
            denom3 <- W_RNA[[x]] %*% (diag(seq(J)) * diag(Pi_RNA[[x]] * diag(H_RNA %*% t(H_RNA))))
        }else{
            denom3 <- 0
        }
        W_RNA[[x]] * ((numer1 + numer2) / (denom1 + denom2 + denom3))^.rho(Beta, root)
    })
}
