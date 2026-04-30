# Low-rank T update: T[k] = U[k] * t(V[k])
# MU updates for U and V derived from:
#   dL/dT = S_hat^(b-1) * G * H * Wt - (S_hat^(b-2) * X) * G * H * Wt
#   dL/dU = dL/dT * V,  dL/dV = t(dL/dT) * U

.updateT2_lowrank <- function(W_RNA, X_Epi, H_Sym, U, V, Beta, L1_T, L2_T, root){
    if(is.matrix(X_Epi)){
        uv <- .updateT2_lowrank_Matrix(W_RNA, X_Epi, H_Sym,
            U, V, Beta, L1_T, L2_T, root)
    }else{
        uv <- .updateT2_lowrank_List(W_RNA, X_Epi, H_Sym,
            U, V, Beta, L1_T, L2_T, root)
    }
    uv
}

.updateT2_lowrank_Matrix <- function(W_RNA, X_Epi, H_Sym,
    U, V, Beta, L1_T, L2_T, root){
    TT <- U %*% t(V)
    G <- TT %*% W_RNA
    S_hat <- G %*% H_Sym %*% t(G)
    # Shared intermediates (avoid materializing l x n gradient)
    WtV <- t(W_RNA) %*% V                              # J x r
    HWtV <- H_Sym %*% WtV                              # J x r
    GHWtV <- G %*% HWtV                                # l x r
    # U update
    numer_U <- (S_hat^(Beta - 2) * X_Epi) %*% GHWtV
    denom_U <- S_hat^(Beta - 1) %*% GHWtV + L1_T + L2_T * U
    U <- U * (numer_U / denom_U)^.rho(Beta, root)
    # Recompute T and G after U update
    TT <- U %*% t(V)
    G <- TT %*% W_RNA
    S_hat <- G %*% H_Sym %*% t(G)
    # V update
    UtSX <- t(U) %*% (S_hat^(Beta - 2) * X_Epi)       # r x l
    GtUtSX <- t(G) %*% t(UtSX)                         # J x r... no
    # Direct: numer_V = t(dL_neg/dT) %*% U = W * H * Gt * (S^(b-2)*X) * U
    SbX_U <- (S_hat^(Beta - 2) * X_Epi) %*% U          # l x r
    GtSbXU <- t(G) %*% SbX_U                            # J x r
    numer_V <- W_RNA %*% H_Sym %*% GtSbXU               # n x r
    Sb_U <- S_hat^(Beta - 1) %*% U                      # l x r
    GtSbU <- t(G) %*% Sb_U                              # J x r
    denom_V <- W_RNA %*% H_Sym %*% GtSbU + L1_T + L2_T * V  # n x r
    V <- V * (numer_V / denom_V)^.rho(Beta, root)
    list(U = U, V = V)
}

.updateT2_lowrank_List <- function(W_RNA, X_Epi, H_Sym,
    U, V, Beta, L1_T, L2_T, root){
    K <- length(X_Epi)
    for(k in seq_len(K)){
        TT_k <- U[[k]] %*% t(V[[k]])
        G <- TT_k %*% W_RNA[[k]]
        S_hat <- G %*% H_Sym %*% t(G)
        # U update
        WtV <- t(W_RNA[[k]]) %*% V[[k]]
        HWtV <- H_Sym %*% WtV
        GHWtV <- G %*% HWtV
        numer_U <- (S_hat^(Beta - 2) * X_Epi[[k]]) %*% GHWtV
        denom_U <- S_hat^(Beta - 1) %*% GHWtV + L1_T + L2_T * U[[k]]
        U[[k]] <- U[[k]] * (numer_U / denom_U)^.rho(Beta, root)
        # Recompute after U update
        TT_k <- U[[k]] %*% t(V[[k]])
        G <- TT_k %*% W_RNA[[k]]
        S_hat <- G %*% H_Sym %*% t(G)
        # V update
        SbX_U <- (S_hat^(Beta - 2) * X_Epi[[k]]) %*% U[[k]]
        GtSbXU <- t(G) %*% SbX_U
        numer_V <- W_RNA[[k]] %*% H_Sym %*% GtSbXU
        Sb_U <- S_hat^(Beta - 1) %*% U[[k]]
        GtSbU <- t(G) %*% Sb_U
        denom_V <- W_RNA[[k]] %*% H_Sym %*% GtSbU + L1_T + L2_T * V[[k]]
        V[[k]] <- V[[k]] * (numer_V / denom_V)^.rho(Beta, root)
    }
    list(U = U, V = V)
}
