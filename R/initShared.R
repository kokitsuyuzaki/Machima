# Initialize w_0 and delta for shared background model

.initShared <- function(X_RNA, X_Epi, W_RNA, T, J, init_g0, init_delta){
    if(is.matrix(X_RNA)){
        .initShared_Matrix(X_RNA, X_Epi, W_RNA, T, J, init_g0, init_delta)
    }else{
        .initShared_List(X_RNA, X_Epi, W_RNA, T, J, init_g0, init_delta)
    }
}

.initShared_Matrix <- function(X_RNA, X_Epi, W_RNA, T, J, init_g0, init_delta){
    n <- nrow(X_RNA)
    # w_0: shared background basis
    if(!is.null(init_g0)){
        w_0 <- init_g0
    }else{
        # Leading singular vector of X_Epi, mapped back through T
        sv <- svd(X_Epi, nu = 1, nv = 1)
        g_0_init <- abs(sv$u[, 1]) * sqrt(sv$d[1])
        # Solve for w_0: g_0 = T * w_0 → w_0 ≈ pinv(T) * g_0
        # Use non-negative least squares approximation: abs(t(T) %*% g_0) normalized
        w_0 <- pmax(as.numeric(t(T) %*% g_0_init), 1e-5)
    }
    # delta: cell-type-specific deviations from W_RNA differential
    if(!is.null(init_delta)){
        delta <- init_delta
    }else{
        W_mean <- rowMeans(W_RNA)
        delta <- pmax(W_RNA - W_mean, 1e-5)
    }
    # h_vec: length J+1, h_0 and h_1..h_J
    h_vec <- rep(1, J + 1)
    list(w_0 = w_0, delta = delta, delta_init = delta, h_vec = h_vec)
}

.initShared_List <- function(X_RNA, X_Epi, W_RNA, T, J, init_g0, init_delta){
    K <- length(X_RNA)
    # w_0: per-chrom shared background basis
    if(!is.null(init_g0)){
        w_0 <- init_g0
    }else{
        w_0 <- lapply(seq_len(K), function(k){
            sv <- svd(X_Epi[[k]], nu = 1, nv = 1)
            g_0_init <- abs(sv$u[, 1]) * sqrt(sv$d[1])
            pmax(as.numeric(t(T[[k]]) %*% g_0_init), 1e-5)
        })
    }
    # delta: per-chrom cell-type-specific deviations
    if(!is.null(init_delta)){
        delta <- init_delta
    }else{
        delta <- lapply(seq_len(K), function(k){
            W_mean <- rowMeans(W_RNA[[k]])
            pmax(W_RNA[[k]] - W_mean, 1e-5)
        })
    }
    # h_vec: shared across chroms
    h_vec <- rep(1, J + 1)
    list(w_0 = w_0, delta = delta, delta_init = delta, h_vec = h_vec)
}
