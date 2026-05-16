.checkMachima2 <- function(X_RNA, X_Epi, label, T,
    fixW_RNA, fixH_RNA, fixT, fixH_Sym,
    orthW_RNA, lambda_orthW, orthH_RNA, orthT, orthH_Sym,
    pseudocount,
    L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
    L1_T, L2_T, L1_H_Sym, L2_H_Sym, orderReg, horizontal,
    J, Beta, root, thr, viz, figdir, num.iter, verbose,
    init_W_RNA, init_H_RNA, init_H_Sym,
    nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
    T_regularization, lambda_T, T_rank, H_Sym_structure,
    lambda_balance,
    J_hic_only, W_hic_init, fixW_hic,
    lambda_coupling, init_U, fixU,
    use_shared_background, init_g0, init_delta, lambda_delta, lambda_delta_anchor, fix_g0){
    # Check X_RNA
    check1 <- is.matrix(X_RNA)
    check2 <- is.list(X_RNA)
    if(!check1 && !check2){
        msg <- paste0("Please specify X_RNA as a matrix or ",
            "a list containing multiple matrices")
        stop(msg)
    }
    if(check1){
        stopifnot(!all(X_RNA == 0))
    }
    if(check2){
        lapply(X_RNA, function(x){
            stopifnot(!all(x == 0))
        })
        nc1 <- ncol(X_RNA[[1]])
        ncs <- unlist(lapply(X_RNA, ncol))
        all(ncs == nc1)
    }
    # Check X_Epi (must be square and symmetric)
    check3 <- is.matrix(X_Epi)
    check4 <- is.list(X_Epi)
    if(!check3 && !check4){
        msg <- paste0("Please specify X_Epi as a symmetric matrix or ",
            "a list containing multiple symmetric matrices")
        stop(msg)
    }
    if(check3){
        stopifnot(!all(X_Epi == 0))
        if(nrow(X_Epi) != ncol(X_Epi)){
            stop("X_Epi must be a square matrix (symmetric)")
        }
        if(!isSymmetric(unname(X_Epi))){
            stop("X_Epi must be a symmetric matrix")
        }
    }
    if(check4){
        lapply(X_Epi, function(x){
            stopifnot(!all(x == 0))
            if(nrow(x) != ncol(x)){
                stop("Each X_Epi must be a square matrix (symmetric)")
            }
            if(!isSymmetric(unname(x))){
                stop("Each X_Epi must be a symmetric matrix")
            }
        })
    }
    # Check X_RNA and X_Epi
    if((check2 && !check4) || (!check2 && check4)){
        msg <- paste0("Please specify both X_RNA and X_Epi as lists")
        stop(msg)
    }
    if(check2 && check4){
        stopifnot(length(X_RNA) == length(X_Epi))
    }
    # Check label
    if(!is.null(label)){
        stopifnot(is.vector(label))
        stopifnot(is.character(label))
        if(check2 && check4){
            stopifnot(length(label) == ncol(X_RNA[[1]]))
        }else{
            stopifnot(length(label) == ncol(X_RNA))
        }
    }
    # Check T
    check5 <- is.matrix(T)
    check6 <- is.null(T)
    check7 <- is.list(T)
    if(!check5 && !check6 && !check7){
        msg <- paste0("Please specify T as a matrix, NULL or ",
            "a list containing multiple matrices")
        stop(msg)
    }
    if(check5){
        stopifnot(check1)
        stopifnot(check3)
        stopifnot(identical(dim(T), c(nrow(X_Epi), nrow(X_RNA))))
    }
    if(check7){
        stopifnot(check2)
        stopifnot(check4)
        lapply(seq_along(T), function(x){
            stopifnot(identical(dim(T[[x]]),
                c(nrow(X_Epi[[x]]), nrow(X_RNA[[x]]))))
        })
    }
    # Check fix
    stopifnot(is.logical(fixW_RNA))
    stopifnot(is.logical(fixH_RNA))
    stopifnot(is.logical(fixT))
    stopifnot(is.logical(fixH_Sym))
    # Deprecation warning for learned T without regularization
    if(!fixT && T_regularization == "none" && check6){
        warning(
            "Learning T (fixT = FALSE) without T_regularization is not ",
            "recommended for paired scATAC + Hi-C with shared bin grids. ",
            "Consider fixT = TRUE (default in v1.2.0+) or ",
            "T_regularization = \"frobenius_unit\" / \"low_rank\".")
    }
    # Check Orthogonal
    stopifnot(is.logical(orthW_RNA))
    stopifnot(is.numeric(lambda_orthW))
    stopifnot(length(lambda_orthW) == 1)
    stopifnot(lambda_orthW >= 0)
    stopifnot(lambda_orthW <= 1)
    stopifnot(is.logical(orthH_RNA))
    stopifnot(is.logical(orthT))
    stopifnot(is.logical(orthH_Sym))
    # Check Pseudo-count
    stopifnot(pseudocount >= 0)
    # Check Regularization Parameters
    stopifnot(L1_W_RNA >= 0)
    stopifnot(L2_W_RNA >= 0)
    stopifnot(L1_H_RNA >= 0)
    stopifnot(L2_H_RNA >= 0)
    stopifnot(L1_T >= 0)
    stopifnot(L2_T >= 0)
    stopifnot(L1_H_Sym >= 0)
    stopifnot(L2_H_Sym >= 0)
    stopifnot(is.logical(orderReg))
    # Check horizontal
    stopifnot(is.logical(horizontal))
    if(horizontal){
        stopifnot(!check6)
    }
    # Check J
    if(check1 && check3){
        stopifnot(J <= min(dim(X_RNA), nrow(X_Epi)))
    }
    if(check2 && check4){
        lapply(seq_along(X_RNA), function(x){
            stopifnot(J <= min(dim(X_RNA[[x]]), nrow(X_Epi[[x]])))
        })
    }
    # Check Beta
    stopifnot(is.numeric(Beta))
    # root
    stopifnot(is.logical(root))
    # Check thr
    stopifnot(thr >= 0)
    # viz
    stopifnot(is.logical(viz))
    # Check figdir
    if(!is.character(figdir) && !is.null(figdir)){
        stop("Please specify the figdir as a string or NULL")
    }
    # Check num.iter
    stopifnot(num.iter >= 0)
    # Check verbose
    stopifnot(is.logical(verbose))
    # Check init_W_RNA
    if(!is.null(init_W_RNA)){
        if(check1){
            if(!is.matrix(init_W_RNA)){
                stop("init_W_RNA must be a matrix when X_RNA is a matrix")
            }
            if(nrow(init_W_RNA) != nrow(X_RNA)){
                stop(paste0("init_W_RNA has ", nrow(init_W_RNA),
                    " rows but X_RNA has ", nrow(X_RNA), " rows"))
            }
            if(ncol(init_W_RNA) != J){
                stop(paste0("init_W_RNA has ", ncol(init_W_RNA),
                    " columns but J = ", J))
            }
        }
        if(check2){
            if(!is.list(init_W_RNA)){
                stop("init_W_RNA must be a list when X_RNA is a list")
            }
            if(length(init_W_RNA) != length(X_RNA)){
                stop(paste0("init_W_RNA has length ", length(init_W_RNA),
                    " but X_RNA has length ", length(X_RNA)))
            }
            lapply(seq_along(init_W_RNA), function(x){
                if(!is.matrix(init_W_RNA[[x]])){
                    stop(paste0("init_W_RNA[[", x, "]] must be a matrix"))
                }
                if(nrow(init_W_RNA[[x]]) != nrow(X_RNA[[x]])){
                    stop(paste0("init_W_RNA[[", x, "]] has ", nrow(init_W_RNA[[x]]),
                        " rows but X_RNA[[", x, "]] has ", nrow(X_RNA[[x]]), " rows"))
                }
                if(ncol(init_W_RNA[[x]]) != J){
                    stop(paste0("init_W_RNA[[", x, "]] has ", ncol(init_W_RNA[[x]]),
                        " columns but J = ", J))
                }
            })
        }
    }
    # Check init_H_RNA
    if(!is.null(init_H_RNA)){
        if(!is.matrix(init_H_RNA)){
            stop("init_H_RNA must be a matrix")
        }
        if(check1){
            expected_m <- ncol(X_RNA)
        }else{
            expected_m <- ncol(X_RNA[[1]])
        }
        if(nrow(init_H_RNA) != J){
            stop(paste0("init_H_RNA has ", nrow(init_H_RNA),
                " rows but J = ", J))
        }
        if(ncol(init_H_RNA) != expected_m){
            stop(paste0("init_H_RNA has ", ncol(init_H_RNA),
                " columns but expected ", expected_m))
        }
    }
    # Check init_H_Sym
    if(!is.null(init_H_Sym)){
        if(!is.matrix(init_H_Sym)){
            stop("init_H_Sym must be a matrix")
        }
        if(nrow(init_H_Sym) != J || ncol(init_H_Sym) != J){
            stop(paste0("init_H_Sym must be ", J, " x ", J,
                " but is ", nrow(init_H_Sym), " x ", ncol(init_H_Sym)))
        }
        if(!isSymmetric(unname(init_H_Sym))){
            stop("init_H_Sym must be a symmetric matrix")
        }
    }
    # Check NMF init parameters
    stopifnot(is.numeric(nmf_init_n_restart))
    stopifnot(nmf_init_n_restart >= 1)
    stopifnot(is.numeric(nmf_init_num_iter))
    stopifnot(nmf_init_num_iter >= 1)
    stopifnot(is.character(nmf_init_algorithm))
    # Check T_regularization
    stopifnot(T_regularization %in% c("none", "frobenius_unit", "l2", "low_rank"))
    stopifnot(is.numeric(lambda_T))
    stopifnot(lambda_T >= 0)
    # Check T_rank
    if(T_regularization == "low_rank"){
        if(is.null(T_rank)){
            stop("T_rank is required when T_regularization = 'low_rank'")
        }
        stopifnot(is.numeric(T_rank))
        stopifnot(T_rank >= 1)
        T_rank <- as.integer(T_rank)
        if(fixT){
            stop("fixT = TRUE is incompatible with T_regularization = 'low_rank'")
        }
        if(check1 && check3){
            if(T_rank > min(nrow(X_Epi), nrow(X_RNA))){
                stop(paste0("T_rank = ", T_rank, " exceeds min(l, n) = ",
                    min(nrow(X_Epi), nrow(X_RNA))))
            }
        }
        if(check2 && check4){
            lapply(seq_along(X_RNA), function(x){
                max_r <- min(nrow(X_Epi[[x]]), nrow(X_RNA[[x]]))
                if(T_rank > max_r){
                    stop(paste0("T_rank = ", T_rank,
                        " exceeds min(l, n) = ", max_r,
                        " for element ", x))
                }
            })
        }
    }else{
        if(!is.null(T_rank)){
            warning("T_rank ignored when T_regularization is not 'low_rank'")
        }
    }
    # Check lambda_balance
    stopifnot(is.numeric(lambda_balance))
    stopifnot(length(lambda_balance) == 1)
    stopifnot(lambda_balance >= 0)
    stopifnot(lambda_balance <= 1)
    # Check H_Sym_structure
    stopifnot(H_Sym_structure %in% c("symmetric", "diagonal"))
    # Check J_hic_only
    stopifnot(is.numeric(J_hic_only))
    stopifnot(length(J_hic_only) == 1)
    stopifnot(as.integer(J_hic_only) >= 0)
    # Check W_hic_init
    if(!is.null(W_hic_init)){
        if(check1){
            if(!is.matrix(W_hic_init)) stop("W_hic_init must be a matrix when X_RNA is a matrix")
            if(nrow(W_hic_init) != nrow(X_Epi)) stop("W_hic_init nrow must match nrow(X_Epi)")
            if(ncol(W_hic_init) != J_hic_only) stop(paste0("W_hic_init ncol must be J_hic_only=", J_hic_only))
            if(any(W_hic_init < 0)) stop("W_hic_init must be non-negative")
        }
        if(check2){
            if(!is.list(W_hic_init)) stop("W_hic_init must be a list when X_RNA is a list")
            if(length(W_hic_init) != length(X_Epi)) stop("W_hic_init length must match X_Epi length")
            lapply(seq_along(W_hic_init), function(k){
                if(nrow(W_hic_init[[k]]) != nrow(X_Epi[[k]])) stop(paste0("W_hic_init[[",k,"]] nrow mismatch"))
                if(ncol(W_hic_init[[k]]) != J_hic_only) stop(paste0("W_hic_init[[",k,"]] ncol must be ",J_hic_only))
                if(any(W_hic_init[[k]] < 0)) stop(paste0("W_hic_init[[",k,"]] must be non-negative"))
            })
        }
    }
    # Check fixW_hic
    stopifnot(is.logical(fixW_hic))
    if(fixW_hic && is.null(W_hic_init) && J_hic_only > 0L){
        stop("fixW_hic=TRUE requires W_hic_init when J_hic_only > 0")
    }
    # Check lambda_coupling
    stopifnot(is.numeric(lambda_coupling))
    stopifnot(length(lambda_coupling) == 1)
    stopifnot(is.infinite(lambda_coupling) || lambda_coupling >= 0)
    # Check init_U
    if(!is.null(init_U)){
        if(check1){
            if(!is.matrix(init_U)) stop("init_U must be a matrix when X_RNA is a matrix")
            if(nrow(init_U) != nrow(X_RNA)) stop("init_U nrow must match nrow(X_RNA)")
            if(ncol(init_U) != J) stop(paste0("init_U ncol must be J=", J))
            if(any(init_U < 0)) stop("init_U must be non-negative")
        }
        if(check2){
            if(!is.list(init_U)) stop("init_U must be a list when X_RNA is a list")
            if(length(init_U) != length(X_RNA)) stop("init_U length must match X_RNA length")
            lapply(seq_along(init_U), function(k){
                if(nrow(init_U[[k]]) != nrow(X_RNA[[k]])) stop(paste0("init_U[[",k,"]] nrow mismatch"))
                if(ncol(init_U[[k]]) != J) stop(paste0("init_U[[",k,"]] ncol must be J=",J))
                if(any(init_U[[k]] < 0)) stop(paste0("init_U[[",k,"]] must be non-negative"))
            })
        }
    }
    # Check/auto-default fixU
    if(is.null(fixU)) fixU <- is.infinite(lambda_coupling)
    stopifnot(is.logical(fixU))
    stopifnot(length(fixU) == 1)
    # Check shared background args
    stopifnot(is.logical(use_shared_background))
    if(use_shared_background && !is.infinite(lambda_coupling)){
        stop("use_shared_background=TRUE is mutually exclusive with lambda_coupling < Inf")
    }
    stopifnot(is.numeric(lambda_delta))
    stopifnot(length(lambda_delta) == 1)
    stopifnot(lambda_delta >= 0 || is.infinite(lambda_delta))
    stopifnot(is.numeric(lambda_delta_anchor))
    stopifnot(length(lambda_delta_anchor) == 1)
    stopifnot(lambda_delta_anchor >= 0 || is.infinite(lambda_delta_anchor))
    stopifnot(is.logical(fix_g0))
}
