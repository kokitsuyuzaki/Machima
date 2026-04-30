.checkMachima2 <- function(X_RNA, X_Epi, label, T,
    fixW_RNA, fixH_RNA, fixT, fixH_Sym,
    orthW_RNA, orthH_RNA, orthT, orthH_Sym,
    pseudocount,
    L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
    L1_T, L2_T, L1_H_Sym, L2_H_Sym, orderReg, horizontal,
    J, Beta, root, thr, viz, figdir, num.iter, verbose,
    init_W_RNA, init_H_RNA, init_H_Sym,
    nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
    T_regularization, lambda_T, T_rank){
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
    # Check Orthogonal
    stopifnot(is.logical(orthW_RNA))
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
}
