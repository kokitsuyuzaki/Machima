# MU update for delta (cell-type-specific deviations epsilon_c)
# With anchoring penalty: lambda_delta_anchor * ||delta_c - delta_c_init||^2
# numer_c += 2 * lambda_delta_anchor * delta_c_init
# denom_c += 2 * lambda_delta_anchor * delta_c + 2 * lambda_delta * delta_c

.updateDelta <- function(X_Epi, T, w_0, delta, delta_init, h_vec, J, Beta,
    L1, L2, lambda_delta, lambda_delta_anchor, root, Pi_Epi){
    if(is.matrix(X_Epi)){
        delta <- .updateDelta_Matrix(X_Epi, T, w_0, delta, delta_init, h_vec, J, Beta,
            L1, L2, lambda_delta, lambda_delta_anchor, root)
    }else{
        delta <- .updateDelta_List(X_Epi, T, w_0, delta, delta_init, h_vec, J, Beta,
            L1, L2, lambda_delta, lambda_delta_anchor, root, Pi_Epi)
    }
    delta
}

.updateDelta_Matrix <- function(X_Epi, T, w_0, delta, delta_init, h_vec, J, Beta,
    L1, L2, lambda_delta, lambda_delta_anchor, root){
    R_full <- .assembleR_shared(T, w_0, delta, h_vec, J)
    A <- R_full^(Beta - 2) * X_Epi
    B <- R_full^(Beta - 1)
    g_0 <- T %*% w_0
    for(c in seq_len(J)){
        d_c <- T %*% delta[, c]
        gd <- g_0 + d_c
        h_c <- h_vec[c + 1]
        numer <- as.numeric(t(T) %*% A %*% (h_c * gd)) +
            2 * lambda_delta_anchor * delta_init[, c]
        denom <- as.numeric(t(T) %*% B %*% (h_c * gd)) +
            L1 + L2 * delta[, c] +
            2 * lambda_delta * delta[, c] +
            2 * lambda_delta_anchor * delta[, c]
        delta[, c] <- delta[, c] * (numer / denom)^.rho(Beta, root)
    }
    delta
}

.updateDelta_List <- function(X_Epi, T, w_0, delta, delta_init, h_vec, J, Beta,
    L1, L2, lambda_delta, lambda_delta_anchor, root, Pi_Epi){
    lapply(seq_along(X_Epi), function(k){
        R_full <- .assembleR_shared(T[[k]], w_0[[k]], delta[[k]], h_vec, J)
        A <- R_full^(Beta - 2) * X_Epi[[k]]
        B <- R_full^(Beta - 1)
        g_0 <- T[[k]] %*% w_0[[k]]
        dk <- delta[[k]]
        dk_init <- delta_init[[k]]
        for(c in seq_len(J)){
            d_c <- T[[k]] %*% dk[, c]
            gd <- g_0 + d_c
            h_c <- h_vec[c + 1]
            numer <- as.numeric(Pi_Epi[[k]] * (t(T[[k]]) %*% A %*% (h_c * gd))) +
                2 * lambda_delta_anchor * dk_init[, c]
            denom <- as.numeric(Pi_Epi[[k]] * (t(T[[k]]) %*% B %*% (h_c * gd))) +
                L1 + L2 * dk[, c] +
                2 * lambda_delta * dk[, c] +
                2 * lambda_delta_anchor * dk[, c]
            dk[, c] <- dk[, c] * (numer / denom)^.rho(Beta, root)
        }
        dk
    })
}
