#' Symmetric Cross-Omics Cell-Type Deconvolution by Joint Non-Negative Matrix Trifactorization
#'
#' For symmetric epigenomic matrices (e.g., similarity/correlation matrices).
#' X_RNA (n x m) is decomposed as W_RNA (n x J) * H_RNA (J x m).
#' X_Epi (l x l, symmetric) is decomposed as G * H_Sym * t(G),
#' where G = T * W_RNA (l x J) and H_Sym is J x J symmetric non-negative.
#'
#' @param X_RNA Single-cell RNA-Seq matrix (n x m) or list of matrices
#' @param X_Epi Symmetric epigenome matrix (l x l) or list of symmetric matrices
#' @param label A length-m character vector to specify the cell type within X_RNA (Default: NULL)
#' @param T Coefficient matrix to connect the dimension of X_RNA and X_Epi (l x n, Default: NULL)
#' @param fixW_RNA Fix value option of W_RNA (Default: FALSE)
#' @param fixH_RNA Fix value option of H_RNA (Default: FALSE)
#' @param fixT Fix value option of T (Default: FALSE)
#' @param fixH_Sym Fix value option of H_Sym (Default: FALSE)
#' @param orthW_RNA Orthogonal option of W_RNA (Default: FALSE)
#' @param orthH_RNA Orthogonal option of H_RNA (Default: FALSE)
#' @param orthT Orthogonal option of T (Default: FALSE)
#' @param orthH_Sym Orthogonal option of H_Sym (Default: FALSE)
#' @param pseudocount Pseudo count to avoid zero division (Default: Machine Epsilon)
#' @param L1_W_RNA Parameter for L1-norm regularization of W_RNA (Default: 1e-10)
#' @param L2_W_RNA Parameter for L2-norm regularization of W_RNA (Default: 1e-10)
#' @param L1_H_RNA Parameter for L1-norm regularization of H_RNA (Default: 1e-10)
#' @param L2_H_RNA Parameter for L2-norm regularization of H_RNA (Default: 1e-10)
#' @param L1_T Parameter for L1-norm regularization of T (Default: 1e-10)
#' @param L2_T Parameter for L2-norm regularization of T (Default: 1e-10)
#' @param L1_H_Sym Parameter for L1-norm regularization of H_Sym (Default: 1e-10)
#' @param L2_H_Sym Parameter for L2-norm regularization of H_Sym (Default: 1e-10)
#' @param orderReg Order regularization (Default: FALSE)
#' @param horizontal Horizontal mode for joint NMF (Default: FALSE)
#' @param J Rank parameter to decompose (Default: 3)
#' @param Beta Parameter of Beta-divergence (Default: 2)
#' @param root Option to add sqrt to the update equation (Default: FALSE)
#' @param thr The threshold to stop the iteration (Default: 1e-10)
#' @param viz Whether the temporal result is visualized (Default: FALSE)
#' @param figdir The figure directory (Default: NULL)
#' @param init Initial value algorithm (Default: "Random")
#' @param num.iter The number of iteration (Default: 30)
#' @param verbose Verbose option (Default: FALSE)
#' @param init_W_RNA User-supplied initial W_RNA: a matrix (n x J) or list of matrices (Default: NULL)
#' @param init_H_RNA User-supplied initial H_RNA: a (J x m) matrix (Default: NULL)
#' @param init_H_Sym User-supplied initial H_Sym: a (J x J) symmetric matrix (Default: NULL)
#' @param nmf_init_n_restart Number of NMF restarts for init (Default: 1)
#' @param nmf_init_num_iter Number of NMF iterations for init (Default: 30)
#' @param nmf_init_algorithm NMF algorithm for init (Default: "Frobenius")
#' @param T_regularization Regularization strategy for T to prevent H_Sym collapse: "none", "frobenius_unit" (per-iteration normalization), or "l2" (penalty). Ignored when fixT=TRUE. (Default: "none")
#' @param lambda_T L2 penalty strength for T when T_regularization="l2" (Default: 0)
#' @return A list containing W_RNA, H_RNA, H_Sym, T, RecError, RelChange
#' @examples
#' X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
#' S <- matrix(runif(15*15), nrow=15, ncol=15)
#' X_Epi <- (S + t(S)) / 2
#' out <- Machima2(X_RNA, X_Epi, T=NULL, verbose=TRUE)
#' \dontrun{
#' # Two-stage workflow: pre-compute NMF, then freeze RNA factors
#' nmf_res <- nnTensor::NMF(X_RNA, J=3, num.iter=200, algorithm="KL")
#' out2 <- Machima2(X_RNA, X_Epi,
#'     init_W_RNA=nmf_res$U, init_H_RNA=t(nmf_res$V),
#'     fixW_RNA=TRUE, fixH_RNA=TRUE, J=3)
#' }
#' @export
Machima2 <- function(X_RNA, X_Epi, label=NULL, T=NULL,
    fixW_RNA=FALSE, fixH_RNA=FALSE, fixT=FALSE, fixH_Sym=FALSE,
    orthW_RNA=FALSE, orthH_RNA=FALSE, orthT=FALSE, orthH_Sym=FALSE,
    pseudocount=.Machine$double.eps,
    L1_W_RNA=1e-10, L2_W_RNA=1e-10,
    L1_H_RNA=1e-10, L2_H_RNA=1e-10,
    L1_T=1e-10, L2_T=1e-10,
    L1_H_Sym=1e-10, L2_H_Sym=1e-10,
    orderReg=FALSE, horizontal=FALSE,
    J=3, Beta=2, root=FALSE, thr=1e-10, viz=FALSE, figdir=NULL,
    init = c("Random", "RandomEpi", "RandomRNA"),
    num.iter=30, verbose=FALSE,
    init_W_RNA=NULL, init_H_RNA=NULL, init_H_Sym=NULL,
    nmf_init_n_restart=1L, nmf_init_num_iter=30L,
    nmf_init_algorithm="Frobenius",
    T_regularization=c("none", "frobenius_unit", "l2"),
    lambda_T=0){
    # Argument Check
    init <- match.arg(init)
    T_regularization <- match.arg(T_regularization)
    .checkMachima2(X_RNA, X_Epi, label, T,
        fixW_RNA, fixH_RNA, fixT, fixH_Sym,
        orthW_RNA, orthH_RNA, orthT, orthH_Sym,
        pseudocount,
        L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
        L1_T, L2_T, L1_H_Sym, L2_H_Sym, orderReg, horizontal,
        J, Beta, root, thr, viz, figdir, num.iter, verbose,
        init_W_RNA, init_H_RNA, init_H_Sym,
        nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
        T_regularization, lambda_T)
    # Initialization
    int <- .initMachima2(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
        init_W_RNA, init_H_RNA, init_H_Sym,
        nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm)
    X_RNA <- int$X_RNA
    X_Epi <- int$X_Epi
    W_RNA <- int$W_RNA
    H_RNA <- int$H_RNA
    H_Sym <- int$H_Sym
    T <- int$T
    Pi_RNA <- int$Pi_RNA
    Pi_Epi <- int$Pi_Epi
    RecError <- int$RecError
    RelChange <- int$RelChange
    # Before Update
    if(viz && !is.null(figdir)){
        png(filename = paste0(figdir, "/0.png"),
            width=2000, height=500)
        .multiImagePlots_Sym(X_RNA, W_RNA, H_RNA, X_Epi, H_Sym, T)
        dev.off()
    }
    if(viz && is.null(figdir)){
        .multiImagePlots_Sym(X_RNA, W_RNA, H_RNA, X_Epi, H_Sym, T)
    }
    # Horizontal Mode: precompute X_GAM (fixed throughout iteration)
    if(horizontal){
        X_GAM <- .updateGAM2_HZL(X_Epi, T)
    }
    # Iteration
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)){
        if(horizontal){
            pre_Error <- .recErrors2_HZL(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi)
        }else{
            pre_Error <- .recErrors2(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi)
        }
        # Horizontal Mode
        if(horizontal){
            # Update1: H_Sym
            if(!fixH_Sym){
                H_Sym <- .updateH_Sym_HZL(X_GAM, W_RNA, H_Sym, J, Beta,
                    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
            }
            # Update2: W_RNA
            if(!fixW_RNA){
                W_RNA <- .updateW_RNA2_HZL(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
                    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
            }
            # Update3: H_RNA
            if(!fixH_RNA){
                H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, J, Beta,
                    L1_H_RNA, L2_H_RNA, orderReg, orthH_RNA, root, Pi_RNA, Pi_Epi)
            }
        # Tri-factorization Mode
        }else{
            # Step1: Update H_Sym
            if(!fixH_Sym){
                H_Sym <- .updateH_Sym(X_Epi, W_RNA, H_Sym, T, J, Beta,
                    L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
            }
            # Step2: Update T
            if(!fixT){
                effective_L2_T <- L2_T
                if(T_regularization == "l2") effective_L2_T <- effective_L2_T + lambda_T
                T <- .updateT2(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, effective_L2_T, orthT, root)
                # Frobenius unit normalization: fix T scale, absorb magnitude into H_Sym
                if(T_regularization == "frobenius_unit"){
                    frob <- .frobNormT(T)
                    T <- .rescaleT(T, frob)
                    H_Sym <- H_Sym * frob$mean_norm_sq
                }
            }
            # Step3: Update W_RNA
            if(!fixW_RNA){
                W_RNA <- .updateW_RNA2(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
                    L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
            }
            # Step4: Update H_RNA
            if(!fixH_RNA){
                H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, J, Beta,
                    L1_H_RNA, L2_H_RNA, orderReg, orthH_RNA, root, Pi_RNA, Pi_Epi)
            }
        }
        # After Update
        if(verbose){
            cat(paste0(iter, " / ", num.iter, "\n"))
        }
        iter <- iter + 1
        if(horizontal){
            RecError[iter] <- .recErrors2_HZL(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi)
        }else{
            RecError[iter] <- .recErrors2(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi)
        }
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if(viz && !is.null(figdir)){
            png(filename = paste0(figdir, "/", iter, ".png"),
                width=2000, height=500)
            .multiImagePlots_Sym(X_RNA, W_RNA, H_RNA, X_Epi, H_Sym, T)
            dev.off()
        }
        if(viz && is.null(figdir)){
            .multiImagePlots_Sym(X_RNA, W_RNA, H_RNA, X_Epi, H_Sym, T)
        }
    }
    # Label Transfer
    if(!is.null(label)){
        asn <- .assignCelltypeNames2(X_RNA, X_Epi, label, W_RNA, H_RNA, H_Sym)
        W_RNA <- asn$W_RNA
        H_RNA <- asn$H_RNA
        H_Sym <- asn$H_Sym
    }
    # Output
    list(W_RNA=W_RNA, H_RNA=H_RNA, H_Sym=H_Sym,
        T=T, RecError=RecError, RelChange=RelChange)
}
