# --- T Frobenius normalization helpers ---
# Uses a single global scalar (RMS of per-chrom norms) in list mode
# to preserve reconstruction exactly for all chromosomes.

.frobNormT <- function(T){
    eps <- sqrt(.Machine$double.eps)
    if(is.matrix(T)){
        scalar <- max(norm(T, "F"), eps)
        list(scalar = scalar, scalar_sq = scalar^2)
    }else{
        norms <- sapply(T, function(t) norm(t, "F"))
        scalar <- max(sqrt(mean(norms^2)), eps)
        list(scalar = scalar, scalar_sq = scalar^2)
    }
}

.rescaleT <- function(T, frob){
    if(is.matrix(T)){
        T / frob$scalar
    }else{
        lapply(T, function(t) t / frob$scalar)
    }
}

# --- Reconstruction errors for symmetric model ---

.recErrors2 <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        d_Beta <- .recErrors2_Matrix(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi)
    }else{
        d_Beta <- .recErrors2_List(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi)
    }
    d_Beta
}

.recErrors2_Matrix <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi){
    G <- T %*% W_RNA
    left <- Pi_RNA * .BetaDivergence(X_RNA, W_RNA %*% H_RNA, Beta)
    right <- Pi_Epi * .BetaDivergence(X_Epi, G %*% H_Sym %*% t(G), Beta)
    left + right
}

.recErrors2_List <- function(X_RNA, W_RNA, H_RNA, X_Epi, T, H_Sym, Beta, Pi_RNA, Pi_Epi){
    lefts <- sum(unlist(lapply(seq_along(X_RNA), function(x){
        Pi_RNA[[x]] * .BetaDivergence(X_RNA[[x]], W_RNA[[x]] %*% H_RNA, Beta)
    })))
    rights <- sum(unlist(lapply(seq_along(X_Epi), function(x){
        G <- T[[x]] %*% W_RNA[[x]]
        Pi_Epi[[x]] * .BetaDivergence(X_Epi[[x]], G %*% H_Sym %*% t(G), Beta)
    })))
    lefts + rights
}

# --- Reconstruction errors for horizontal mode ---

.recErrors2_HZL <- function(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi){
    if(is.matrix(X_RNA)){
        d_Beta <- .recErrors2_HZL_Matrix(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi)
    }else{
        d_Beta <- .recErrors2_HZL_List(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi)
    }
    d_Beta
}

.recErrors2_HZL_Matrix <- function(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi){
    left <- Pi_RNA * .BetaDivergence(X_RNA, W_RNA %*% H_RNA, Beta)
    right <- Pi_Epi * .BetaDivergence(X_GAM, W_RNA %*% H_Sym %*% t(W_RNA), Beta)
    left + right
}

.recErrors2_HZL_List <- function(X_RNA, W_RNA, H_RNA, X_GAM, H_Sym, Beta, Pi_RNA, Pi_Epi){
    lefts <- sum(unlist(lapply(seq_along(X_RNA), function(x){
        Pi_RNA[[x]] * .BetaDivergence(X_RNA[[x]], W_RNA[[x]] %*% H_RNA, Beta)
    })))
    rights <- sum(unlist(lapply(seq_along(X_GAM), function(x){
        Pi_Epi[[x]] * .BetaDivergence(X_GAM[[x]], W_RNA[[x]] %*% H_Sym %*% t(W_RNA[[x]]), Beta)
    })))
    lefts + rights
}

# --- X_GAM for horizontal mode ---

.updateGAM2_HZL <- function(X_Epi, T){
    if(is.matrix(X_Epi)){
        t(T) %*% X_Epi %*% T
    }else{
        lapply(seq_along(X_Epi), function(x){
            t(T[[x]]) %*% X_Epi[[x]] %*% T[[x]]
        })
    }
}

# --- Cell type name assignment for symmetric model ---
# Uses shared .estimateCelltypes() from assignCelltypeNames.R

.assignCelltypeNames2 <- function(X_RNA, X_Epi, label, W_RNA, H_RNA, H_Sym){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        asn <- .assignCelltypeNames2_Matrix(X_RNA, label, W_RNA, H_RNA, H_Sym)
    }else{
        asn <- .assignCelltypeNames2_List(X_RNA, label, W_RNA, H_RNA, H_Sym)
    }
    asn
}

.assignCelltypeNames2_Matrix <- function(X_RNA, label, W_RNA, H_RNA, H_Sym){
    estimated.celltypes <- .estimateCelltypes(X_RNA, label, W_RNA)
    colnames(W_RNA) <- estimated.celltypes
    rownames(H_RNA) <- estimated.celltypes
    rownames(H_Sym) <- estimated.celltypes
    colnames(H_Sym) <- estimated.celltypes
    list(W_RNA=W_RNA, H_RNA=H_RNA, H_Sym=H_Sym)
}

.assignCelltypeNames2_List <- function(X_RNA, label, W_RNA, H_RNA, H_Sym){
    X_RNA2 <- do.call("rbind", X_RNA)
    W_RNA2 <- do.call("rbind", W_RNA)
    asn <- .assignCelltypeNames2_Matrix(X_RNA2, label, W_RNA2, H_RNA, H_Sym)
    list(W_RNA=.Mat2List(X_RNA, asn$W_RNA),
        H_RNA=asn$H_RNA, H_Sym=asn$H_Sym)
}

# --- Visualization for symmetric model ---

.multiImagePlots_Sym <- function(X_RNA, W_RNA, H_RNA, X_Epi, H_Sym, T){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        G <- T %*% W_RNA
        recX_RNA <- W_RNA %*% H_RNA
        recX_Epi <- G %*% H_Sym %*% t(G)
        layout(rbind(1:4, 5:8))
        image.plot2(X_RNA, main="X_RNA")
        image.plot2(recX_RNA, main="rec X_RNA")
        image.plot2(W_RNA, main="W_RNA")
        image.plot2(H_RNA, main="H_RNA")
        image.plot2(X_Epi, main="X_Epi")
        image.plot2(recX_Epi, main="rec X_Epi")
        image.plot2(H_Sym, main="H_Sym")
        image.plot2(T, main="T")
    }else{
        X_RNA2 <- do.call("rbind", X_RNA)
        W_RNA2 <- do.call("rbind", W_RNA)
        recX_RNA2 <- W_RNA2 %*% H_RNA
        layout(rbind(1:4, 5:8))
        image.plot2(X_RNA2, main="X_RNA")
        image.plot2(recX_RNA2, main="rec X_RNA")
        image.plot2(W_RNA2, main="W_RNA")
        image.plot2(H_RNA, main="H_RNA")
        image.plot2(do.call("rbind", X_Epi), main="X_Epi")
        recX_Epi <- do.call("rbind",
            lapply(seq_along(T), function(x){
            G <- T[[x]] %*% W_RNA[[x]]
            G %*% H_Sym %*% t(G)
        }))
        image.plot2(recX_Epi, main="rec X_Epi")
        image.plot2(H_Sym, main="H_Sym")
        plot.new()
    }
}
