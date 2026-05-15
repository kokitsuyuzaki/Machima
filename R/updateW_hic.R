# MU update for W_hic (Hi-C-only basis columns)
# W_hic[[k]] lives in Hi-C bin space (l_k × J_hic_only)
# R_full = G_shared·H_Sym·G_sharedᵀ + W_hic·diag(h_hic)·W_hicᵀ

.updateW_hic <- function(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic, Beta,
    L1_W, L2_W, root, Pi_Epi, U=NULL){
    if(is.matrix(X_Epi)){
        W_hic <- .updateW_hic_Matrix(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
            Beta, L1_W, L2_W, root, U_k=U)
    }else{
        W_hic <- .updateW_hic_List(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
            Beta, L1_W, L2_W, root, Pi_Epi, U=U)
    }
    W_hic
}

.updateW_hic_Matrix <- function(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
    Beta, L1_W, L2_W, root, U_k=NULL){
    R_full <- .reconstructEpi_single(W_RNA, T, W_hic, H_Sym, h_hic, U_k=U_k)
    numer <- (R_full^(Beta - 2) * X_Epi) %*% W_hic %*% diag(h_hic, nrow=length(h_hic))
    denom <- R_full^(Beta - 1) %*% W_hic %*% diag(h_hic, nrow=length(h_hic)) + L1_W + L2_W * W_hic
    W_hic * (numer / denom)^.rho(Beta, root)
}

.updateW_hic_List <- function(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
    Beta, L1_W, L2_W, root, Pi_Epi, U=NULL){
    lapply(seq_along(X_Epi), function(k){
        U_k <- if(!is.null(U)) U[[k]] else NULL
        R_full <- .reconstructEpi_single(W_RNA[[k]], T[[k]], W_hic[[k]], H_Sym, h_hic, U_k=U_k)
        Dh <- diag(h_hic, nrow=length(h_hic))
        numer <- Pi_Epi[[k]] * ((R_full^(Beta - 2) * X_Epi[[k]]) %*% W_hic[[k]] %*% Dh)
        denom <- Pi_Epi[[k]] * (R_full^(Beta - 1) %*% W_hic[[k]] %*% Dh + L1_W + L2_W * W_hic[[k]])
        W_hic[[k]] * (numer / denom)^.rho(Beta, root)
    })
}
