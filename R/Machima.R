#' Cross-Omics Cell-Type Deconvolution by Joint Non-Negative Matrix Trifactorization
#' Gene expression matrix and epigenomic feature matrix are decomposed at once
#' by non-negative matrix trifactorization algorithm.
#' For the details of the method, see the reference section of GitHub README.md
#' <https://github.com/kokitsuyuzaki/Machima>.
#' @param X_RNA Single-cell RNA-Seq matrix (n x m)
#' @param X_Epi Bulk Epigenome feature matrix (l x o)
#' @param label A length-m character vector to specify the cell type within X_RNA (Default: NULL)
#' @param T Coefficient matrix to connect the dimension of X_RNA and X_Epi (l x n, Default: NULL)
#' @param fixW_RNA Fix value option of W_RNA (for Transfer Learning, Default: FALSE)
#' @param fixH_RNA Fix value option of H_RNA (for Transfer Learning, Default: FALSE)
#' @param fixT Fix value option of T (for Transfer Learning, Default: FALSE)
#' @param orthW_RNA Orthogonal option of W_RNA (for uniqueness, Default: FALSE)
#' @param orthH_RNA Orthogonal option of H_RNA (for uniqueness, Default: FALSE)
#' @param orthT Orthogonal option of T (for uniqueness, Default: FALSE)
#' @param orthH_Epi Orthogonal option of H_Epi (for uniqueness, Default: FALSE)
#' @param pseudocount Pseudo count to avoid zero division (Default: Machine Epsilon)
#' @param L1_W_RNA Parameter for L1-norm regularization of W_RNA (Default: 1e-10)
#' @param L2_W_RNA Parameter for L2-norm regularization of W_RNA (Default: 1e-10)
#' @param L1_H_RNA Parameter for L1-norm regularization of H_RNA (Default: 1e-10)
#' @param L2_H_RNA Parameter for L2-norm regularization of H_RNA (Default: 1e-10)
#' @param L1_T Parameter for L1-norm regularization of T (Default: 1e-10)
#' @param L2_T Parameter for L2-norm regularization of T (Default: 1e-10)
#' @param L1_H_Epi Parameter for L1-norm regularization of H_Epi (Default: 1e-10)
#' @param L2_H_Epi Parameter for L2-norm regularization of H_Epi (Default: 1e-10)
#' @param orderReg Order regularization to sort the column vectors of W_RNA, H_RNA, and H_Epi by the L2 norm in ascending order (Default: FALSE)
#' @param J Rank parameter to decompose (Default: 3)
#' @param Beta Parameter of Beta-divergence (Default: 2)
#' @param root Option to add sqrt to the update equation (Default: FALSE)
#' @param thr The threshold to stop the iteration (Default: 1e-10)
#' @param viz Whether the temporal result is visualized (Default: FALSE)
#' @param figdir The figure directory (Default: NULL)
#' @param init Initial value algorithm (Default: "RandomEpi")
#' @param num.iter The number of iteration (Default: 30)
#' @param verbose Verbose option (Default: FALSE)
#' @return A list containing the Joint-nonnegative Trifactorization
#' @examples
#' X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
#' X_Epi <- matrix(runif(15*25), nrow=15, ncol=25)
#' out <- Machima(X_RNA, X_Epi, T=NULL, verbose=TRUE)
#' @importFrom nnTensor NMF
#' @importFrom fields image.plot
#' @importFrom graphics layout
#' @importFrom grDevices dev.off png
#' @importFrom stats runif cor
#' @importFrom igraph graph_from_incidence_matrix max_bipartite_match
#' @export
Machima <- function(X_RNA, X_Epi, label=NULL, T=NULL,
    fixW_RNA=FALSE, fixH_RNA=FALSE, fixT=FALSE,
    orthW_RNA=FALSE, orthH_RNA=FALSE, orthT=FALSE, orthH_Epi=FALSE,
    pseudocount=.Machine$double.eps,
    L1_W_RNA=1e-10, L2_W_RNA=1e-10,
    L1_H_RNA=1e-10, L2_H_RNA=1e-10,
    L1_T=1e-10, L2_T=1e-10,
    L1_H_Epi=1e-10, L2_H_Epi=1e-10,
    orderReg=FALSE,
    J=3, Beta=2, root=FALSE, thr=1e-10, viz=FALSE, figdir=NULL,
    init = c("RandomEpi", "RandomRNA", "Random", "NMFAlign", "NMFAlign2"),
    num.iter=30, verbose=FALSE){
    # Argument Check
    init <- match.arg(init)
    .checkMachima(X_RNA, X_Epi, label, T,
        fixW_RNA, fixH_RNA, fixT,
        orthW_RNA, orthH_RNA, orthT, orthH_Epi,
        pseudocount,
        L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
        L1_T, L2_T, L1_H_Epi, L2_H_Epi, orderReg,
        J, Beta, root, thr, viz, figdir, num.iter, verbose)
    # Initialization
    int <- .initMachima(X_RNA, X_Epi, T, fixT, pseudocount, J, init, thr)
    X_RNA <- int$X_RNA
    X_Epi <- int$X_Epi
    W_RNA <- int$W_RNA
    H_RNA <- int$H_RNA
    H_Epi <- int$H_Epi
    T <- int$T
    Pi_RNA <- int$Pi_RNA
    Pi_Epi <- int$Pi_Epi
    RecError <- int$RecError
    RelChange <- int$RelChange
    X_GAM <- .updateGAM(X_RNA, X_Epi, W_RNA, H_Epi)
    # Before Update
    if(viz && !is.null(figdir)){
        png(filename = paste0(figdir, "/0.png"),
            width=2500, height=500)
        .multiImagePlots3(X_RNA, W_RNA, H_RNA, X_GAM, X_Epi, H_Epi, T)
        dev.off()
    }
    if(viz && is.null(figdir)){
        .multiImagePlots3(X_RNA, W_RNA, H_RNA, X_GAM, X_Epi, H_Epi, T)
    }
    # Iteration
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)){
        pre_Error <- .recErrors(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi)
        # Step1: Update H_Epi
        H_Epi <- .updateH_Epi(X_Epi, W_RNA, H_Epi, T, J, Beta, L1_H_Epi, L2_H_Epi, orderReg, orthH_Epi, root, Pi_RNA, Pi_Epi)
        # Step2: Update T
        if(!fixT){
            T <- .updateT(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T, orthT, root)
        }
        # Step3: Update W_RNA
        if(!fixW_RNA){
            W_RNA <- .updateW_RNA(X_RNA, X_Epi, W_RNA, H_RNA, H_Epi, T, J, Beta, L1_W_RNA, L2_W_RNA, orderReg, orthW_RNA, root, Pi_RNA, Pi_Epi)
        }
        # Step4: Update H_RNA
        if(!fixH_RNA){
            H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, J, Beta,
                L1_H_RNA, L2_H_RNA, orderReg, orthH_RNA, root, Pi_RNA, Pi_Epi)
        }
        # X_GAM for visualization
        X_GAM <- .updateGAM(X_RNA, X_Epi, W_RNA, H_Epi)
        # After Update
        if(verbose){
            cat(paste0(iter, " / ", num.iter, "\n"))
        }
        iter <- iter + 1
        RecError[iter] <- .recErrors(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Epi, Beta, Pi_RNA, Pi_Epi)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if(viz && !is.null(figdir)){
            png(filename = paste0(figdir, "/", iter, ".png"),
                width=2500, height=500)
            .multiImagePlots3(X_RNA, W_RNA, H_RNA, X_GAM, X_Epi, H_Epi, T)
            dev.off()
        }
        if(viz && is.null(figdir)){
            .multiImagePlots3(X_RNA, W_RNA, H_RNA, X_GAM, X_Epi, H_Epi, T)
        }
    }
    # Label Transfer
    if(!is.null(label)){
        asn <- .assignCelltypeNames(X_RNA, X_Epi, label, W_RNA, H_RNA, H_Epi)
        W_RNA <- asn$W_RNA
        H_RNA <- asn$H_RNA
        H_Epi <- asn$H_Epi
    }
    # Output
    list(W_RNA=W_RNA, H_RNA=H_RNA, H_Epi=H_Epi, X_GAM=X_GAM,
        T=T, RecError=RecError, RelChange=RelChange)
}
