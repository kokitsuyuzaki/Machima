# MU update for delta (cell-type-specific deviations epsilon_c)
# num_c = T^T * A_k * h_c * (g_0 + delta_c)
# denom_c = T^T * B_k * h_c * (g_0 + delta_c) + L1 + L2*eps_c + 2*lambda_delta*eps_c

.updateDelta <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1, L2, lambda_delta, root, Pi_Epi){
    if(is.matrix(X_Epi)){
        delta <- .updateDelta_Matrix(X_Epi, T, w_0, delta, h_vec, J, Beta,
            L1, L2, lambda_delta, root)
    }else{
        delta <- .updateDelta_List(X_Epi, T, w_0, delta, h_vec, J, Beta,
            L1, L2, lambda_delta, root, Pi_Epi)
    }
    delta
}

.updateDelta_Matrix <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1, L2, lambda_delta, root){
    R_full <- .assembleR_shared(T, w_0, delta, h_vec, J)
    A <- R_full^(Beta - 2) * X_Epi
    B <- R_full^(Beta - 1)
    g_0 <- T %*% w_0
    for(c in seq_len(J)){
        d_c <- T %*% delta[, c]
        gd <- g_0 + d_c
        h_c <- h_vec[c + 1]
        numer <- t(T) %*% A %*% (h_c * gd)
        denom <- t(T) %*% B %*% (h_c * gd) + L1 + L2 * delta[, c] +
            2 * lambda_delta * delta[, c]
        delta[, c] <- delta[, c] * (as.numeric(numer) / as.numeric(denom))^.rho(Beta, root)
    }
    delta
}

.updateDelta_List <- function(X_Epi, T, w_0, delta, h_vec, J, Beta,
    L1, L2, lambda_delta, root, Pi_Epi){
    lapply(seq_along(X_Epi), function(k){
        R_full <- .assembleR_shared(T[[k]], w_0[[k]], delta[[k]], h_vec, J)
        A <- R_full^(Beta - 2) * X_Epi[[k]]
        B <- R_full^(Beta - 1)
        g_0 <- T[[k]] %*% w_0[[k]]
        dk <- delta[[k]]
        for(c in seq_len(J)){
            d_c <- T[[k]] %*% dk[, c]
            gd <- g_0 + d_c
            h_c <- h_vec[c + 1]
            numer <- Pi_Epi[[k]] * (t(T[[k]]) %*% A %*% (h_c * gd))
            denom <- Pi_Epi[[k]] * (t(T[[k]]) %*% B %*% (h_c * gd)) +
                L1 + L2 * dk[, c] + 2 * lambda_delta * dk[, c]
            dk[, c] <- dk[, c] * (as.numeric(numer) / as.numeric(denom))^.rho(Beta, root)
        }
        dk
    })
}
