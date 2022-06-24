#' Cross-omics and chromosome Cell-type Deconvolution by Joint Non-negative Trifactorization
#' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' @param X_RNA Single-cell RNA-Seq matrix (n x m)
#' @param X_Epi Bulk Epigenome feature matrix (l x o)
#' @param fixW_RNA Fix value option of W_RNA (for Transfer Learning)
#' @param fixH_RNA Fix value option of H_RNA (for Transfer Learning)
#' @param fixT Fix value option of T (for Transfer Learning)
#' @param pseudocount Pseudo count
#' @param L1_W_RNA Parameter for L1-norm regularization of W_RNA
#' @param L2_W_RNA Parameter for L2-norm regularization of W_RNA
#' @param L1_H_RNA Parameter for L1-norm regularization of H_RNA
#' @param L2_H_RNA Parameter for L2-norm regularization of H_RNA
#' @param L1_T Parameter for L1-norm regularization of T
#' @param L2_T Parameter for L2-norm regularization of T
#' @param L1_H_Epi Parameter for L1-norm regularization of H_Epi
#' @param L2_H_Epi Parameter for L2-norm regularization of H_Epi
#' @param J Rank parameter to decompose
#' @param Beta Parameter of Beta-divergence
#' @param thr The threshold to stop the iteration
#' @param viz Whether the temporal result is visualized (Default: FALSE)
#' @param figdir The figure directory
#' @param init Initial value algorithm
#' @param num.iter The number of iteration
#' @param verbose Verbose option
#' @return
#' @examples
#' X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
#' X_Epi <- matrix(runif(15*25), nrow=15, ncol=25)
#' out <- Machima(X_RNA, X_Epi, T=NULL, verbose=TRUE)
#' @importFrom nnTensor NMF
#' @importFrom fields image.plot
#' @export
Machima <- function(X_RNA, X_Epi, T=NULL,
    fixW_RNA=TRUE, fixH_RNA=TRUE, fixT=FALSE,
    pseudocount=1e-10,
    L1_W_RNA=1e-10, L2_W_RNA=1e-10,
    L1_H_RNA=1e-10, L2_H_RNA=1e-10,
    L1_T=1e-10, L2_T=1e-10,
    L1_H_Epi=1e-10, L2_H_Epi=1e-10,
    J=3, Beta=2, thr=1e-10, viz=FALSE, figdir=NULL,
    init = c("NMF", "Random"),
    num.iter=30, verbose=FALSE){
    # Argument Check
    init <- match.arg(init)
    .checkJointBetaNMTF(X_RNA, X_Epi, T,
        fixW_RNA, fixH_RNA, fixT, pseudocount,
        L1_W_RNA, L2_W_RNA, L1_H_RNA, L2_H_RNA,
        L1_T, L2_T, L1_H_Epi, L2_H_Epi,
        J, Beta, thr, viz, figdir, num.iter, verbose)
    # Initialization
    int <- .initJointBetaNMTF(X_RNA, X_Epi, T, fixT, pseudocount, J, init)
    X_RNA <- int$X_RNA
    X_Epi <- int$X_Epi
    W_RNA <- int$W_RNA
    H_RNA <- int$H_RNA
    W_Epi <- int$W_Epi
    H_Epi <- int$H_Epi
    T <- int$T
    X_GAM <- .updateGAM(X_RNA, X_Epi, W_RNA, H_Epi)
    # Before Update
    if(viz && !is.null(figdir)){
        png(filename = paste0(figdir, "/0.png"),
            width=2500, height=500)
        .multiImagePlots3(X_RNA, X_Epi, X_GAM)
        dev.off()
    }
    if(viz && is.null(figdir)){
        .multiImagePlots3(X_RNA, X_Epi, X_GAM)
    }
    # Iteration
    iter <- 1
    while (iter <= num.iter){
        # Step1: Update W_RNA
        if(!fixW_RNA){
            W_RNA <- .updateW_RNA(X_RNA, W_RNA, H_RNA, H_Epi, T, Beta, L1_W_RNA, L2_W_RNA)
        }
        # Step2: Update H_RNA
        if(!fixH_RNA){
            H_RNA <- .updateH_RNA(X_RNA, W_RNA, H_RNA, Beta, L1_H_RNA, L2_H_RNA)
        }
        # Step3: Update T
        if(!fixT){
            T <- .updateT(W_RNA, X_Epi, H_Epi, T, Beta, L1_T, L2_T)
        }
        # Step4: Update H_Epi
        H_Epi <- .updateH_Epi(X_Epi, W_RNA, H_Epi, T, Beta, L1_H_Epi, L2_H_Epi)
        # X_GAM
        X_GAM <- .updateGAM(X_RNA, X_Epi, W_RNA, H_Epi)
        # After Update
        if(viz && !is.null(figdir)){
            png(filename = paste0(figdir, "/", iter, ".png"),
                width=2500, height=500)
            .multiImagePlots3(X_RNA, X_Epi, X_GAM)
            dev.off()
        }
        if(viz && is.null(figdir)){
            .multiImagePlots3(X_RNA, X_Epi, X_GAM)
        }
        if(verbose){
            cat(paste0(iter, " / ", num.iter, "\n"))
        }
        iter <- iter + 1
    }
    # Output
    list(W_RNA=W_RNA, H_RNA=H_RNA,
        W_Epi=W_Epi, H_Epi=H_Epi,
        X_GAM=X_GAM, T=T)
}
