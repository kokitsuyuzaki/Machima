# MU update for U (Hi-C-specific deviation from W_RNA)
# W_E = W_RNA + U; Hi-C loss gradient + coupling penalty 2*lambda*U in denom

.updateU <- function(X_Epi, W_RNA, U, H_Sym, T, Beta,
    L1_W_RNA, L2_W_RNA, lambda_coupling, root, Pi_Epi,
    W_hic=NULL, h_hic=numeric(0)){
    if(is.matrix(X_Epi)){
        U <- .updateU_Matrix(X_Epi, W_RNA, U, H_Sym, T, Beta,
            L1_W_RNA, L2_W_RNA, lambda_coupling, root, Pi_Epi,
            W_hic, h_hic)
    }else{
        U <- .updateU_List(X_Epi, W_RNA, U, H_Sym, T, Beta,
            L1_W_RNA, L2_W_RNA, lambda_coupling, root, Pi_Epi,
            W_hic, h_hic)
    }
    U
}

.updateU_Matrix <- function(X_Epi, W_RNA, U, H_Sym, T, Beta,
    L1_W_RNA, L2_W_RNA, lambda_coupling, root, Pi_Epi,
    W_hic, h_hic){
    W_E <- W_RNA + U
    G <- T %*% W_E
    GH <- G %*% H_Sym
    S_hat <- GH %*% t(G)
    S_hat <- .addHic(S_hat, W_hic, h_hic)
    numer <- Pi_Epi * (t(T) %*% (S_hat^(Beta - 2) * X_Epi) %*% GH)
    denom <- Pi_Epi * (t(T) %*% S_hat^(Beta - 1) %*% GH) +
        L1_W_RNA + L2_W_RNA * U + 2 * lambda_coupling * U
    U * (numer / denom)^.rho(Beta, root)
}

.updateU_List <- function(X_Epi, W_RNA, U, H_Sym, T, Beta,
    L1_W_RNA, L2_W_RNA, lambda_coupling, root, Pi_Epi,
    W_hic, h_hic){
    lapply(seq_along(X_Epi), function(k){
        W_E <- W_RNA[[k]] + U[[k]]
        G <- T[[k]] %*% W_E
        GH <- G %*% H_Sym
        S_hat <- GH %*% t(G)
        W_hic_k <- if(!is.null(W_hic)) W_hic[[k]] else NULL
        S_hat <- .addHic(S_hat, W_hic_k, h_hic)
        numer <- Pi_Epi[[k]] * (t(T[[k]]) %*% (S_hat^(Beta - 2) * X_Epi[[k]]) %*% GH)
        denom <- Pi_Epi[[k]] * (t(T[[k]]) %*% S_hat^(Beta - 1) %*% GH) +
            L1_W_RNA + L2_W_RNA * U[[k]] + 2 * lambda_coupling * U[[k]]
        U[[k]] * (numer / denom)^.rho(Beta, root)
    })
}
