# MU update for w_0 (shared background basis)
# num = T^T * A_k * (h_0*g_0 + sum_c h_c*(g_0+delta_c))
# denom = T^T * B_k * (...same...) + L1 + L2*w_0

.updateSharedG <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1, L2, root, Pi_Epi){
    if(is.matrix(X_Epi)){
        w_0 <- .updateSharedG_Matrix(X_Epi, T, w_0, delta, h_vec, J, Beta,
            L1, L2, root)
    }else{
        w_0 <- .updateSharedG_List(X_Epi, T, w_0, delta, h_vec, J, Beta,
            L1, L2, root, Pi_Epi)
    }
    w_0
}

.updateSharedG_Matrix <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1, L2, root){
    R_full <- .assembleR_shared(T, w_0, delta, h_vec, J)
    A <- R_full^(Beta - 2) * X_Epi
    B <- R_full^(Beta - 1)
    # Weighted sum: h_0*g_0 + sum_c h_c*(g_0+delta_c)
    g_0 <- T %*% w_0
    wsum <- h_vec[1] * g_0
    for(c in seq_len(J)){
        d_c <- T %*% delta[, c]
        wsum <- wsum + h_vec[c + 1] * (g_0 + d_c)
    }
    numer <- t(T) %*% A %*% wsum
    denom <- t(T) %*% B %*% wsum + L1 + L2 * w_0
    w_0 * (as.numeric(numer) / as.numeric(denom))^.rho(Beta, root)
}

.updateSharedG_List <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1, L2, root, Pi_Epi){
    lapply(seq_along(X_Epi), function(k){
        R_full <- .assembleR_shared(T[[k]], w_0[[k]], delta[[k]], h_vec, J)
        A <- R_full^(Beta - 2) * X_Epi[[k]]
        B <- R_full^(Beta - 1)
        g_0 <- T[[k]] %*% w_0[[k]]
        wsum <- h_vec[1] * g_0
        for(c in seq_len(J)){
            d_c <- T[[k]] %*% delta[[k]][, c]
            wsum <- wsum + h_vec[c + 1] * (g_0 + d_c)
        }
        numer <- Pi_Epi[[k]] * (t(T[[k]]) %*% A %*% wsum)
        denom <- Pi_Epi[[k]] * (t(T[[k]]) %*% B %*% wsum) + L1 + L2 * w_0[[k]]
        w_0[[k]] * (as.numeric(numer) / as.numeric(denom))^.rho(Beta, root)
    })
}
