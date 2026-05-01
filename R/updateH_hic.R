# Diagonal MU update for h_hic (Hi-C-only diagonal weights)
# Same form as updateH_Sym_diag but using W_hic instead of G_shared

.updateH_hic <- function(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic, Beta,
    L1_H, L2_H, root, Pi_Epi){
    if(is.matrix(X_Epi)){
        h_hic <- .updateH_hic_Matrix(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
            Beta, L1_H, L2_H, root)
    }else{
        h_hic <- .updateH_hic_List(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
            Beta, L1_H, L2_H, root, Pi_Epi)
    }
    h_hic
}

.updateH_hic_Matrix <- function(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
    Beta, L1_H, L2_H, root){
    Jh <- length(h_hic)
    R_full <- .reconstructEpi_single(W_RNA, T, W_hic, H_Sym, h_hic)
    SbX <- R_full^(Beta - 2) * X_Epi
    Sb <- R_full^(Beta - 1)
    numer <- rep(0, Jh)
    denom <- rep(0, Jh)
    for(i in seq_len(Jh)){
        wi <- W_hic[, i]
        numer[i] <- as.numeric(t(wi) %*% SbX %*% wi)
        denom[i] <- as.numeric(t(wi) %*% Sb %*% wi) + L1_H + L2_H * h_hic[i]
    }
    h_hic * (numer / denom)^.rho(Beta, root)
}

.updateH_hic_List <- function(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic,
    Beta, L1_H, L2_H, root, Pi_Epi){
    Jh <- length(h_hic)
    numer <- rep(0, Jh)
    denom <- rep(0, Jh)
    for(k in seq_along(X_Epi)){
        R_full <- .reconstructEpi_single(W_RNA[[k]], T[[k]], W_hic[[k]], H_Sym, h_hic)
        SbX <- R_full^(Beta - 2) * X_Epi[[k]]
        Sb <- R_full^(Beta - 1)
        for(i in seq_len(Jh)){
            wi <- W_hic[[k]][, i]
            numer[i] <- numer[i] + Pi_Epi[[k]] * as.numeric(t(wi) %*% SbX %*% wi)
            denom[i] <- denom[i] + Pi_Epi[[k]] * (as.numeric(t(wi) %*% Sb %*% wi) + L1_H + L2_H * h_hic[i])
        }
    }
    h_hic * (numer / denom)^.rho(Beta, root)
}
