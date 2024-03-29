.checkMachima <- function(X_RNA, X_Epi, label, T,
    fixW_RNA, fixH_RNA, fixT,
    orthW_RNA, orthH_RNA, orthT, orthH_Epi,
    pseudocount,
    L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
    L1_T, L2_T, L1_H_Epi, L2_H_Epi, orderReg, horizontal,
    J, Beta, root, thr, viz, figdir, num.iter, verbose){
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
    # Check X_Epi
    check3 <- is.matrix(X_Epi)
    check4 <- is.list(X_Epi)
    if(!check3 && !check4){
        msg <- paste0("Please specify X_Epi as a matrix or ",
            "a list containing multiple matrices")
        stop(msg)
    }
    if(check3){
        stopifnot(!all(X_Epi == 0))
    }
    if(check4){
        lapply(X_Epi, function(x){
            stopifnot(!all(x == 0))
        })
        nc1 <- ncol(X_Epi[[1]])
        ncs <- unlist(lapply(X_Epi, ncol))
        all(ncs == nc1)
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
    # Check Orthogonal
    stopifnot(is.logical(orthW_RNA))
    stopifnot(is.logical(orthH_RNA))
    stopifnot(is.logical(orthT))
    stopifnot(is.logical(orthH_Epi))
    # Check Pseudo-count
    stopifnot(pseudocount >= 0)
    # Check Regularization Parameters
    stopifnot(L1_W_RNA >= 0)
    stopifnot(L2_W_RNA >= 0)
    stopifnot(L1_H_RNA >= 0)
    stopifnot(L2_H_RNA >= 0)
    stopifnot(L1_T >= 0)
    stopifnot(L2_T >= 0)
    stopifnot(L1_H_Epi >= 0)
    stopifnot(L2_H_Epi >= 0)
    stopifnot(is.logical(orderReg))
    # Check horizontal
    stopifnot(is.logical(horizontal))
    if(horizontal){
        stopifnot(!check6)
    }
    # Check J
    if(check1 && check3){
        stopifnot(J <= min(dim(X_RNA), dim(X_Epi)))
    }
    if(check2 && check4){
        lapply(seq_along(X_RNA), function(x){
            stopifnot(J <= min(dim(X_RNA[[x]]), dim(X_Epi[[x]])))
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
}
