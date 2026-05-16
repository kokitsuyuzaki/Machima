# MU update for h_vec (J+1 diagonal weights in shared background model)
# h_vec[1] = h_0 (background), h_vec[2..J+1] = h_1..h_J (cell-type)
# For each i: basis_i is g_0 (if i=1) or g_0+delta_{i-1} (if i>1)
# num_i = sum_k basis_i^T * A_k * basis_i
# denom_i = sum_k basis_i^T * B_k * basis_i + L1 + L2*h_i

.updateH_shared <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1_H, L2_H, root, Pi_Epi){
    if(is.matrix(X_Epi)){
        h_vec <- .updateH_shared_Matrix(X_Epi, T, w_0, delta, h_vec, J, Beta,
            L1_H, L2_H, root)
    }else{
        h_vec <- .updateH_shared_List(X_Epi, T, w_0, delta, h_vec, J, Beta,
            L1_H, L2_H, root, Pi_Epi)
    }
    h_vec
}

.updateH_shared_Matrix <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1_H, L2_H, root){
    Jp1 <- J + 1
    R_full <- .assembleR_shared(T, w_0, delta, h_vec, J)
    A <- R_full^(Beta - 2) * X_Epi
    B <- R_full^(Beta - 1)
    g_0 <- T %*% w_0
    numer <- rep(0, Jp1)
    denom <- rep(0, Jp1)
    # h_0 (background): basis = g_0
    numer[1] <- as.numeric(t(g_0) %*% A %*% g_0)
    denom[1] <- as.numeric(t(g_0) %*% B %*% g_0) + L1_H + L2_H * h_vec[1]
    # h_c (cell-type): basis = g_0 + delta_c
    for(c in seq_len(J)){
        d_c <- T %*% delta[, c]
        gd <- g_0 + d_c
        numer[c + 1] <- as.numeric(t(gd) %*% A %*% gd)
        denom[c + 1] <- as.numeric(t(gd) %*% B %*% gd) + L1_H + L2_H * h_vec[c + 1]
    }
    h_vec * (numer / denom)^.rho(Beta, root)
}

.updateH_shared_List <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1_H, L2_H, root, Pi_Epi){
    Jp1 <- J + 1
    numer <- rep(0, Jp1)
    denom <- rep(0, Jp1)
    for(k in seq_along(X_Epi)){
        R_full <- .assembleR_shared(T[[k]], w_0[[k]], delta[[k]], h_vec, J)
        A <- R_full^(Beta - 2) * X_Epi[[k]]
        B <- R_full^(Beta - 1)
        g_0 <- T[[k]] %*% w_0[[k]]
        numer[1] <- numer[1] + Pi_Epi[[k]] * as.numeric(t(g_0) %*% A %*% g_0)
        denom[1] <- denom[1] + Pi_Epi[[k]] * (as.numeric(t(g_0) %*% B %*% g_0) + L1_H + L2_H * h_vec[1])
        for(c in seq_len(J)){
            d_c <- T[[k]] %*% delta[[k]][, c]
            gd <- g_0 + d_c
            numer[c + 1] <- numer[c + 1] + Pi_Epi[[k]] * as.numeric(t(gd) %*% A %*% gd)
            denom[c + 1] <- denom[c + 1] + Pi_Epi[[k]] * (as.numeric(t(gd) %*% B %*% gd) + L1_H + L2_H * h_vec[c + 1])
        }
    }
    h_vec * (numer / denom)^.rho(Beta, root)
}
