#' Symmetric Cross-Omics Cell-Type Deconvolution by Joint Non-Negative Matrix Trifactorization
#'
#' For symmetric epigenomic matrices (e.g., similarity/correlation matrices).
#' X_RNA (n x m) is decomposed as W_RNA (n x J) * H_RNA (J x m).
#' X_Epi (l x l, symmetric) is decomposed as G * H_Sym * t(G),
#' where G = T * W_RNA (l x J) and H_Sym is J x J symmetric non-negative.
#'
#' @details
#' When \code{label} is supplied, dimnames of \code{H_Sym} and rownames of
#' \code{H_RNA} are assigned by correlating each component's H_RNA row with
#' per-cluster average expression, then solving a max-weight bipartite matching
#' (via \code{igraph::max_bipartite_match}). Components that remain unmatched
#' (e.g. when J > number of clusters) are labeled \code{"comp_j"} instead of NA.
#' Such components often act as background-absorbing factors that capture the
#' bulk variance of X_Epi; their H_Sym diagonal entry is typically the largest.
#' Downstream analyses may drop or baseline-correct \code{"comp_*"} entries.
#'
#' When \code{H_Sym_structure = "diagonal"}, the model becomes symmetric CP:
#' \code{X_Epi[k] = sum_i h_i * g_i[k] * g_i[k]^T}, where each component has
#' an independent contact pattern weighted by \code{h_i = H_Sym[i,i]}. This is
#' the natural choice when cross-type "interaction" (off-diagonal) is not
#' biologically meaningful. The default \code{"symmetric"} mode allows off-diagonal
#' entries to capture co-occurrence of contacts between components.
#'
#' For paired scATAC + Hi-C on the same genomic bin grid, T is naturally
#' the identity matrix (the two modalities share the same coordinate axis).
#' Learning T as a free matrix introduces a fictitious coordinate
#' transformation that degrades both reconstruction quality and cell-type
#' discrimination. The v1.2.0 default is therefore \code{fixT = TRUE}.
#'
#' @param X_RNA Single-cell RNA-Seq matrix (n x m) or list of matrices
#' @param X_Epi Symmetric epigenome matrix (l x l) or list of symmetric matrices
#' @param label A length-m character vector to specify the cell type within X_RNA (Default: NULL)
#' @param T Coefficient matrix to connect the dimension of X_RNA and X_Epi (l x n, Default: NULL)
#' @param fixW_RNA Fix value option of W_RNA (Default: FALSE)
#' @param fixH_RNA Fix value option of H_RNA (Default: FALSE)
#' @param fixT If TRUE (default), T is fixed during iteration; when T is also NULL, an identity matrix is auto-constructed per chrom. If FALSE, T is learned as a free dense matrix. For paired scATAC + Hi-C on a shared bin grid, fixT=TRUE is recommended. (Default: TRUE)
#' @param fixH_Sym Fix value option of H_Sym (Default: FALSE)
#' @param orthW_RNA Deprecated. Use lambda_orthW instead. If TRUE, sets lambda_orthW=1. (Default: FALSE)
#' @param lambda_orthW Strength of W_RNA column orthogonality: 0=standard NMF, 1=full orthogonal, intermediate=blend. (Default: 0)
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
#' @param lambda_balance Balance between X_RNA and X_Epi loss terms: 0=RNA only, 1=Epi only, 0.5=equal (default; matches pre-1.3.0). (Default: 0.5)
#' @param T_regularization Regularization strategy for T: "none", "frobenius_unit", "l2", or "low_rank". Ignored when fixT=TRUE. (Default: "none")
#' @param lambda_T L2 penalty strength for T when T_regularization="l2" (Default: 0)
#' @param T_rank Rank of low-rank T parametrization (T=U*t(V)). Required when T_regularization="low_rank". (Default: NULL)
#' @param lambda_coupling Coupling strength between W_RNA and Hi-C basis. Inf=hard share (default), 0=independent. When finite, W_E=W_RNA+U is used for Hi-C with penalty lambda_coupling*||U||^2. (Default: Inf)
#' @param init_U Optional list of initial U matrices (each n_k x J, non-negative). When NULL, U is initialized to small random values runif(1e-5, 1e-2). (Default: NULL)
#' @param fixU If TRUE, U is not updated. Auto-defaults to TRUE when lambda_coupling=Inf. (Default: NULL = auto)
#' @param use_shared_background If TRUE, decompose Hi-C into shared background g_0 + cell-type deviations delta_c. Mutually exclusive with lambda_coupling < Inf. (Default: FALSE)
#' @param init_g0 Optional list of initial w_0 vectors (each length n_k). NULL = SVD-derived. (Default: NULL)
#' @param init_delta Optional list of initial delta matrices (each n_k x J). NULL = W_RNA differential. (Default: NULL)
#' @param lambda_delta Sparsity penalty on delta toward zero. (Default: 0)
#' @param lambda_delta_anchor Anchoring penalty pulling delta toward its RNA-derived init. 0=no anchor, Inf=freeze at init. (Default: 1)
#' @param fix_g0 If TRUE, w_0 is not updated. (Default: FALSE)
#' @param J_hic_only Number of Hi-C-only basis columns. When >0, Hi-C reconstruction becomes G_full*H_full*G_full^T with hic-only columns independent of W_RNA. (Default: 0)
#' @param W_hic_init Optional list of initial W_hic matrices (each l_k x J_hic_only). (Default: NULL)
#' @param fixW_hic If TRUE, do not update W_hic. (Default: FALSE)
#' @param H_Sym_structure Constrain H_Sym to "symmetric" (default; full J(J+1)/2 parameters; symmetric Tucker) or "diagonal" (J diagonal entries only; symmetric CP / sym-PARAFAC). The diagonal option is recommended for cell-type deconvolution of bulk Hi-C. (Default: "symmetric")
#' @return A list containing W_RNA, H_RNA, H_Sym, T, RecError, RelChange. When T_regularization="low_rank", also T_factors. When J_hic_only>0, also W_hic and h_hic.
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
    fixW_RNA=FALSE, fixH_RNA=FALSE, fixT=TRUE, fixH_Sym=FALSE,
    orthW_RNA=FALSE, lambda_orthW=0,
    orthH_RNA=FALSE, orthT=FALSE, orthH_Sym=FALSE,
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
    lambda_balance=0.5,
    T_regularization=c("none", "frobenius_unit", "l2", "low_rank"),
    lambda_T=0, T_rank=NULL,
    H_Sym_structure=c("symmetric", "diagonal"),
    J_hic_only=0L, W_hic_init=NULL, fixW_hic=FALSE,
    lambda_coupling=Inf, init_U=NULL, fixU=NULL,
    use_shared_background=FALSE, init_g0=NULL, init_delta=NULL,
    lambda_delta=0, lambda_delta_anchor=1, fix_g0=FALSE){
    # Argument Check
    init <- match.arg(init)
    T_regularization <- match.arg(T_regularization)
    H_Sym_structure <- match.arg(H_Sym_structure)
    # Auto-default fixU
    if(is.null(fixU)) fixU <- is.infinite(lambda_coupling)
    # Deprecation: orthW_RNA -> lambda_orthW
    if(orthW_RNA){
        warning("orthW_RNA=TRUE is deprecated; setting lambda_orthW=1. Use lambda_orthW directly.")
        lambda_orthW <- 1
    }
    .checkMachima2(X_RNA, X_Epi, label, T,
        fixW_RNA, fixH_RNA, fixT, fixH_Sym,
        orthW_RNA, lambda_orthW, orthH_RNA, orthT, orthH_Sym,
        pseudocount,
        L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
        L1_T, L2_T, L1_H_Sym, L2_H_Sym, orderReg, horizontal,
        J, Beta, root, thr, viz, figdir, num.iter, verbose,
        init_W_RNA, init_H_RNA, init_H_Sym,
        nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
        T_regularization, lambda_T, T_rank, H_Sym_structure,
        lambda_balance, J_hic_only, W_hic_init, fixW_hic,
        lambda_coupling, init_U, fixU,
        use_shared_background, init_g0, init_delta, lambda_delta, lambda_delta_anchor, fix_g0)
    # Initialization
    int <- .initMachima2(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr,
        init_W_RNA, init_H_RNA, init_H_Sym,
        nmf_init_n_restart, nmf_init_num_iter, nmf_init_algorithm,
        T_regularization, T_rank, H_Sym_structure, lambda_balance,
        J_hic_only, W_hic_init, lambda_coupling, init_U)
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
    if(T_regularization == "low_rank"){
        U <- int$U
        V <- int$V
    }
    W_hic <- int$W_hic
    h_hic <- int$h_hic
    U_coupling <- int$U_coupling
    # Shared background init
    shared_w0 <- NULL
    shared_delta <- NULL
    shared_h_vec <- NULL
    if(use_shared_background){
        sh <- .initShared(X_RNA, X_Epi, W_RNA, T, J, init_g0, init_delta)
        shared_w0 <- sh$w_0
        shared_delta <- sh$delta
        shared_delta_init <- sh$delta_init
        shared_h_vec <- sh$h_vec
    }
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
        if(use_shared_background){
            pre_Error <- .recErrors_shared(X_RNA, W_RNA, H_RNA, X_Epi, T,
                shared_w0, shared_delta, shared_h_vec, J, Beta, Pi_RNA, Pi_Epi)
        }else if(horizontal){
            pre_Error <- .recErrors2_HZL(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi)
        }else{
            pre_Error <- .recErrors2(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi, W_hic, h_hic, U_coupling)
        }
        # Shared Background Mode
        if(use_shared_background){
            # Update h_vec (J+1 diagonal weights)
            shared_h_vec <- .updateH_shared(X_Epi, T, shared_w0, shared_delta,
                shared_h_vec, J, Beta, L1_H_Sym, L2_H_Sym, root, Pi_Epi)
            # Update w_0 (shared background basis)
            if(!fix_g0){
                shared_w0 <- .updateSharedG(X_Epi, T, shared_w0, shared_delta,
                    shared_h_vec, J, Beta, L1_W_RNA, L2_W_RNA, root, Pi_Epi)
            }
            # Update delta (cell-type deviations)
            if(!is.infinite(lambda_delta_anchor) && !is.infinite(lambda_delta)){
                # Update delta (skip if either penalty is Inf)
                shared_delta <- .updateDelta(X_Epi, T, shared_w0, shared_delta,
                    shared_delta_init, shared_h_vec, J, Beta,
                    L1_W_RNA, L2_W_RNA, lambda_delta, lambda_delta_anchor,
                    root, Pi_Epi)
            }
            # Update W_RNA and H_RNA (RNA side only, no Hi-C gradient)
            if(!fixW_RNA){
                # RNA-only W_RNA update (standard NMF, no Epi term)
                if(is.matrix(X_RNA)){
                    WH <- W_RNA %*% H_RNA
                    numer <- (WH^(Beta - 2) * X_RNA) %*% t(H_RNA)
                    denom <- WH^(Beta - 1) %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA
                    W_RNA <- W_RNA * (numer / denom)^.rho(Beta, root)
                }else{
                    W_RNA <- lapply(seq_along(X_RNA), function(k){
                        WH <- W_RNA[[k]] %*% H_RNA
                        numer <- Pi_RNA[[k]] * ((WH^(Beta - 2) * X_RNA[[k]]) %*% t(H_RNA))
                        denom <- Pi_RNA[[k]] * (WH^(Beta - 1) %*% t(H_RNA) + L1_W_RNA + L2_W_RNA * W_RNA[[k]])
                        W_RNA[[k]] * (numer / denom)^.rho(Beta, root)
                    })
                }
            }
            if(!fixH_RNA){
                H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, J, Beta,
                    L1_H_RNA, L2_H_RNA, orderReg, orthH_RNA, root, Pi_RNA, Pi_Epi)
            }
        # Horizontal Mode
        }else if(horizontal){
            # Update1: H_Sym
            if(!fixH_Sym){
                if(H_Sym_structure == "diagonal"){
                    h <- .updateH_Sym_diag_HZL(X_GAM, W_RNA, diag(H_Sym), J, Beta,
                        L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi)
                    dn <- dimnames(H_Sym)
                    H_Sym <- diag(h, nrow=J)
                    dimnames(H_Sym) <- dn
                }else{
                    H_Sym <- .updateH_Sym_HZL(X_GAM, W_RNA, H_Sym, J, Beta,
                        L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi)
                }
            }
            # Update2: W_RNA
            if(!fixW_RNA){
                W_RNA <- .updateW_RNA2_HZL(X_RNA, X_GAM, W_RNA, H_RNA, H_Sym, J, Beta,
                    L1_W_RNA, L2_W_RNA, orderReg, lambda_orthW, root, Pi_RNA, Pi_Epi)
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
                if(H_Sym_structure == "diagonal"){
                    h <- .updateH_Sym_diag(X_Epi, W_RNA, diag(H_Sym), T, J, Beta,
                        L1_H_Sym, L2_H_Sym, orderReg, root, Pi_Epi,
                        W_hic=W_hic, h_hic=h_hic, U=U_coupling)
                    dn <- dimnames(H_Sym)
                    H_Sym <- diag(h, nrow=J)
                    dimnames(H_Sym) <- dn
                }else{
                    H_Sym <- .updateH_Sym(X_Epi, W_RNA, H_Sym, T, J, Beta,
                        L1_H_Sym, L2_H_Sym, orderReg, orthH_Sym, root, Pi_RNA, Pi_Epi,
                        W_hic=W_hic, h_hic=h_hic, U=U_coupling)
                }
            }
            # Step2: Update T
            if(!fixT){
                if(T_regularization == "low_rank"){
                    uv <- .updateT2_lowrank(W_RNA, X_Epi, H_Sym, U, V, Beta, L1_T, L2_T, root,
                        W_hic=W_hic, h_hic=h_hic, U_coup=U_coupling)
                    U <- uv$U
                    V <- uv$V
                    if(is.matrix(X_Epi)){
                        T <- U %*% t(V)
                    }else{
                        T <- lapply(seq_along(U), function(k) U[[k]] %*% t(V[[k]]))
                    }
                }else{
                    effective_L2_T <- L2_T
                    if(T_regularization == "l2") effective_L2_T <- effective_L2_T + lambda_T
                    T <- .updateT2(W_RNA, X_Epi, H_Sym, T, Beta, L1_T, effective_L2_T, orthT, root,
                        W_hic=W_hic, h_hic=h_hic, U=U_coupling)
                    if(T_regularization == "frobenius_unit"){
                        frob <- .frobNormT(T)
                        T <- .rescaleT(T, frob)
                        H_Sym <- H_Sym * frob$scalar_sq
                    }
                }
            }
            # Step3: Update W_RNA
            if(!fixW_RNA){
                W_RNA <- .updateW_RNA2(X_RNA, X_Epi, W_RNA, H_RNA, H_Sym, T, J, Beta,
                    L1_W_RNA, L2_W_RNA, orderReg, lambda_orthW, root, Pi_RNA, Pi_Epi,
                    W_hic=W_hic, h_hic=h_hic, U=U_coupling)
            }
            # Step4: Update H_RNA
            if(!fixH_RNA){
                H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, J, Beta,
                    L1_H_RNA, L2_H_RNA, orderReg, orthH_RNA, root, Pi_RNA, Pi_Epi)
            }
            # Step5: Update hic-only factors
            if(J_hic_only > 0L){
                h_hic <- .updateH_hic(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic, Beta,
                    L1_H_Sym, L2_H_Sym, root, Pi_Epi, U=U_coupling)
                if(!fixW_hic){
                    W_hic <- .updateW_hic(X_Epi, W_RNA, W_hic, H_Sym, T, h_hic, Beta,
                        L1_W_RNA, L2_W_RNA, root, Pi_Epi, U=U_coupling)
                }
            }
            # Step6: Update U (soft-coupling deviation)
            if(!is.null(U_coupling) && !fixU){
                U_coupling <- .updateU(X_Epi, W_RNA, U_coupling, H_Sym, T, Beta,
                    L1_W_RNA, L2_W_RNA, lambda_coupling, root, Pi_Epi,
                    W_hic=W_hic, h_hic=h_hic)
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
                    if(use_shared_background){
            RecError[iter] <- .recErrors_shared(X_RNA, W_RNA, H_RNA, X_Epi, T,
                shared_w0, shared_delta, shared_h_vec, J, Beta, Pi_RNA, Pi_Epi)
        }else{
            RecError[iter] <- .recErrors2(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi, W_hic, h_hic, U_coupling)
        }
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
    out <- list(W_RNA=W_RNA, H_RNA=H_RNA, H_Sym=H_Sym,
        T=T, RecError=RecError, RelChange=RelChange)
    if(T_regularization == "low_rank"){
        out$T_factors <- list(U=U, V=V)
    }
    if(J_hic_only > 0L){
        out$W_hic <- W_hic
        out$h_hic <- h_hic
    }
    if(!is.null(U_coupling)){
        out$U <- U_coupling
    }
    if(use_shared_background){
        out$w_0 <- shared_w0
        out$delta <- shared_delta
        out$h_vec <- shared_h_vec
        out$predict_celltype_R <- function(c, chrom_idx=1){
            if(is.matrix(X_Epi)){
                .predict_celltype_R(T, shared_w0, shared_delta, shared_h_vec, c)
            }else{
                .predict_celltype_R(T[[chrom_idx]], shared_w0[[chrom_idx]],
                    shared_delta[[chrom_idx]], shared_h_vec, c)
            }
        }
        out$predict_celltype_delta <- function(c, chrom_idx=1){
            if(is.matrix(X_Epi)){
                .predict_celltype_delta(T, shared_w0, shared_delta, shared_h_vec, c, J)
            }else{
                .predict_celltype_delta(T[[chrom_idx]], shared_w0[[chrom_idx]],
                    shared_delta[[chrom_idx]], shared_h_vec, c, J)
            }
        }
    }
    out
}
