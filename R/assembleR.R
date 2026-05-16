# Assemble R_full for shared background + differential model
# R_full = h_0 * g_0 * g_0^T + sum_c h_c * (g_0 + delta_c) * (g_0 + delta_c)^T

.assembleR_shared <- function(T_k, w_0_k, delta_k, h_vec, J){
    # h_vec: length J+1 where h_vec[1] = h_0, h_vec[2..J+1] = h_1..h_J
    g_0 <- T_k %*% w_0_k
    R <- h_vec[1] * g_0 %*% t(g_0)
    for(c in seq_len(J)){
        d_c <- T_k %*% delta_k[, c]
        gd <- g_0 + d_c
        R <- R + h_vec[c + 1] * gd %*% t(gd)
    }
    R
}

# Per-cell-type prediction: R_c = h_c * (g_0 + delta_c) * (g_0 + delta_c)^T
.predict_celltype_R <- function(T_k, w_0_k, delta_k, h_vec, c){
    g_0 <- T_k %*% w_0_k
    d_c <- T_k %*% delta_k[, c]
    gd <- g_0 + d_c
    h_vec[c + 1] * gd %*% t(gd)
}

# Differential: Delta_c = R_c - R_avg
.predict_celltype_delta <- function(T_k, w_0_k, delta_k, h_vec, c, J){
    R_c <- .predict_celltype_R(T_k, w_0_k, delta_k, h_vec, c)
    h_ct <- h_vec[2:(J + 1)]
    R_avg <- matrix(0, nrow(R_c), ncol(R_c))
    g_0 <- T_k %*% w_0_k
    for(j in seq_len(J)){
        d_j <- T_k %*% delta_k[, j]
        gd <- g_0 + d_j
        R_avg <- R_avg + h_ct[j] * gd %*% t(gd)
    }
    R_avg <- R_avg / sum(h_ct)
    R_c - R_avg
}

# RecError for shared background model
.recErrors_shared <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, w_0, delta,
    h_vec, J, Beta, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        left <- Pi_RNA * .BetaDivergence(X_RNA, W_RNA %*% H_RNA, Beta)
        R_full <- .assembleR_shared(T, w_0, delta, h_vec, J)
        right <- Pi_Epi * .BetaDivergence(X_Epi, R_full, Beta)
        left + right
    }else{
        lefts <- sum(unlist(lapply(seq_along(X_RNA), function(k){
            Pi_RNA[[k]] * .BetaDivergence(X_RNA[[k]], W_RNA[[k]] %*% H_RNA, Beta)
        })))
        rights <- sum(unlist(lapply(seq_along(X_Epi), function(k){
            R_full <- .assembleR_shared(T[[k]], w_0[[k]], delta[[k]], h_vec, J)
            Pi_Epi[[k]] * .BetaDivergence(X_Epi[[k]], R_full, Beta)
        })))
        lefts + rights
    }
}
