.initMachima2 <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
    init_W_RNA, init_H_RNA, init_H_Sym,
    nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
    T_regularization, T_rank){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        int <- .initMachima2_Matrix(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
            init_W_RNA, init_H_RNA, init_H_Sym,
            nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
            T_regularization, T_rank)
    }else{
        int <- .initMachima2_List(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
            init_W_RNA, init_H_RNA, init_H_Sym,
            nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
            T_regularization, T_rank)
    }
    int
}

.initSymH <- function(J){
    A <- matrix(runif(J * J), nrow=J, ncol=J)
    (A + t(A)) / 2
}

.initMachima2_Matrix <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
    init_W_RNA, init_H_RNA, init_H_Sym,
    nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
    T_regularization, T_rank){
    X_RNA[which(X_RNA == 0)] <- pseudocount
    X_Epi[which(X_Epi == 0)] <- pseudocount
    # Symmetrize after pseudocount
    X_Epi <- (X_Epi + t(X_Epi)) / 2
    if(init == "RandomEpi"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA, J=J,
            n_restart=nmf_init_n_restart, num_iter=nmf_init_num_iter,
            algorithm=nmf_init_algorithm), X_RNA)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- t(t(out1_1$U) * out1_2$normA)
        H_RNA <- t(out1_2$A)
        # Random H_Sym (symmetric)
        H_Sym <- .initSymH(J)
        # Random T
        if(!fixT){
            nr <- nrow(X_Epi)
            nc <- nrow(X_RNA)
            T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
        }
    }
    if(init == "RandomRNA"){
        # Random W_RNA/H_RNA
        W_RNA <- .normalizeCols(
            matrix(runif(nrow(X_RNA)*J),
            nrow=nrow(X_RNA), ncol=J))$A
        H_RNA <- t(W_RNA) %*% X_RNA
        # Random H_Sym (symmetric)
        H_Sym <- .initSymH(J)
        # Random T
        if(!fixT){
            nr <- nrow(X_Epi)
            nc <- nrow(X_RNA)
            T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
        }
    }
    if(init == "Random"){
        # Random W_RNA/H_RNA
        W_RNA <- matrix(runif(nrow(X_RNA)*J), nrow=nrow(X_RNA), ncol=J)
        H_RNA <- t(W_RNA) %*% X_RNA
        # Random H_Sym (symmetric)
        H_Sym <- .initSymH(J)
        # Random T
        if(!fixT){
            nr <- nrow(X_Epi)
            nc <- nrow(X_RNA)
            T <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
        }
    }
    # Override with user-supplied init values
    if(!is.null(init_W_RNA)) W_RNA <- init_W_RNA
    if(!is.null(init_H_RNA)) H_RNA <- init_H_RNA
    if(!is.null(init_H_Sym)) H_Sym <- init_H_Sym
    # Low-rank T initialization
    U <- NULL
    V <- NULL
    if(T_regularization == "low_rank"){
        if(is.null(T) || all(T == 0)){
            U <- matrix(runif(nrow(X_Epi) * T_rank), nrow(X_Epi), T_rank)
            V <- matrix(runif(nrow(X_RNA) * T_rank), nrow(X_RNA), T_rank)
        }else{
            sv <- svd(T, nu = T_rank, nv = T_rank)
            sqrtd <- sqrt(sv$d[seq_len(T_rank)])
            U <- abs(sv$u[, seq_len(T_rank), drop=FALSE]) * rep(sqrtd, each=nrow(sv$u))
            V <- abs(sv$v[, seq_len(T_rank), drop=FALSE]) * rep(sqrtd, each=nrow(sv$v))
        }
        T <- U %*% t(V)
    }
    # Weight
    Pi_RNA <- .weight(X_RNA)
    Pi_Epi <- .weight(X_Epi)
    # Error
    RecError <- c()
    RelChange <- c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, H_Sym=H_Sym,
        T=T, U=U, V=V, Pi_RNA=Pi_RNA, Pi_Epi=Pi_Epi,
        RecError=RecError, RelChange=RelChange)
}

.initMachima2_List <- function(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
    init_W_RNA, init_H_RNA, init_H_Sym,
    nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
    T_regularization, T_rank){
    X_RNA <- lapply(X_RNA, function(x){
        x[which(x == 0)] <- pseudocount
        x
    })
    X_Epi <- lapply(X_Epi, function(x){
        x[which(x == 0)] <- pseudocount
        (x + t(x)) / 2
    })
    X_RNA2 <- do.call("rbind", X_RNA)
    if(init == "RandomEpi"){
        # NMF with RNA
        out1_1 <- .reArrangeOuts(.returnBestNMF(X_RNA2, J=J,
            n_restart=nmf_init_n_restart, num_iter=nmf_init_num_iter,
            algorithm=nmf_init_algorithm), X_RNA2)
        out1_2 <- .normalizeCols(out1_1$V)
        W_RNA <- .Mat2List(X_RNA, out1_1$U)
        H_RNA <- t(out1_2$A * out1_2$normA)
        # Random H_Sym (symmetric)
        H_Sym <- .initSymH(J)
        # Random T
        if(!fixT){
            T <- lapply(seq_along(X_Epi), function(x){
                nr <- nrow(X_Epi[[x]])
                nc <- nrow(X_RNA[[x]])
                matrix(runif(nr*nc), nrow=nr, ncol=nc)
            })
        }
    }
    if(init == "RandomRNA"){
        # Random W_RNA/H_RNA
        W_RNA2 <- matrix(runif(nrow(X_RNA2)*J), nrow=nrow(X_RNA2), ncol=J)
        W_RNA <- .Mat2List(X_RNA, W_RNA2)
        H_RNA <- t(W_RNA2) %*% X_RNA2
        # Random H_Sym (symmetric)
        H_Sym <- .initSymH(J)
        # Random T
        if(!fixT){
            T <- lapply(seq_along(X_Epi), function(x){
                nr <- nrow(X_Epi[[x]])
                nc <- nrow(X_RNA[[x]])
                matrix(runif(nr*nc), nrow=nr, ncol=nc)
            })
        }
    }
    if(init == "Random"){
        # Random W_RNA/H_RNA
        W_RNA2 <- matrix(runif(nrow(X_RNA2)*J), nrow=nrow(X_RNA2), ncol=J)
        W_RNA <- .Mat2List(X_RNA, W_RNA2)
        H_RNA <- t(W_RNA2) %*% X_RNA2
        # Random H_Sym (symmetric)
        H_Sym <- .initSymH(J)
        # Random T
        if(!fixT){
            T <- lapply(seq_along(X_Epi), function(x){
                nr <- nrow(X_Epi[[x]])
                nc <- nrow(X_RNA[[x]])
                matrix(runif(nr*nc), nrow=nr, ncol=nc)
            })
        }
    }
    # Override with user-supplied init values
    if(!is.null(init_W_RNA)) W_RNA <- init_W_RNA
    if(!is.null(init_H_RNA)) H_RNA <- init_H_RNA
    if(!is.null(init_H_Sym)) H_Sym <- init_H_Sym
    # Low-rank T initialization
    U <- NULL
    V <- NULL
    if(T_regularization == "low_rank"){
        if(is.null(T) || !is.list(T)){
            U <- lapply(seq_along(X_Epi), function(k){
                matrix(runif(nrow(X_Epi[[k]]) * T_rank), nrow(X_Epi[[k]]), T_rank)
            })
            V <- lapply(seq_along(X_RNA), function(k){
                matrix(runif(nrow(X_RNA[[k]]) * T_rank), nrow(X_RNA[[k]]), T_rank)
            })
        }else{
            U <- lapply(T, function(t){
                sv <- svd(t, nu = T_rank, nv = T_rank)
                sqrtd <- sqrt(sv$d[seq_len(T_rank)])
                abs(sv$u[, seq_len(T_rank), drop=FALSE]) * rep(sqrtd, each=nrow(sv$u))
            })
            V <- lapply(T, function(t){
                sv <- svd(t, nu = T_rank, nv = T_rank)
                sqrtd <- sqrt(sv$d[seq_len(T_rank)])
                abs(sv$v[, seq_len(T_rank), drop=FALSE]) * rep(sqrtd, each=nrow(sv$v))
            })
        }
        T <- lapply(seq_along(U), function(k) U[[k]] %*% t(V[[k]]))
    }
    # Weight
    Pi_RNA <- lapply(X_RNA, .weight)
    Pi_Epi <- lapply(X_Epi, .weight)
    # Error
    RecError <- c()
    RelChange <- c()
    RecError[1] <- thr * 10
    RelChange[1] <- thr * 10
    list(X_RNA=X_RNA, X_Epi=X_Epi,
        W_RNA=W_RNA, H_RNA=H_RNA, H_Sym=H_Sym,
        T=T, U=U, V=V, Pi_RNA=Pi_RNA, Pi_Epi=Pi_Epi,
        RecError=RecError, RelChange=RelChange)
}
