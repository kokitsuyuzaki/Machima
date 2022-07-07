.assignCelltypeNames <- function(X_RNA, X_Epi, label, W_RNA, H_RNA, H_Epi){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        asn <- .assignCelltypeNames_Matrix(X_RNA, label, W_RNA, H_RNA, H_Epi)
    }else{
        asn <- .assignCelltypeNames_List(X_RNA, label, W_RNA, H_RNA, H_Epi)
    }
    asn
}

.assignCelltypeNames_Matrix <- function(X_RNA, label, W_RNA, H_RNA, H_Epi){
    # Setting
    uniq.celltypes <- unique(label)
    components <- paste0("component", seq(ncol(W_RNA)))
    # Average vectors
    AvgX_RNA <- do.call("cbind",
        lapply(uniq.celltypes, function(x){
            target <- which(label == x)
            if(length(target) >= 2){
                rowMeans(X_RNA[, target])
            }else{
                X_RNA[, target]
            }}))
    # Matching by Correlation Coefficient
    cor.matrix <- cor(AvgX_RNA, W_RNA)
    cor.matrix[which(is.na(cor.matrix))] <- 0
    rownames(cor.matrix) <- uniq.celltypes
    colnames(cor.matrix) <- components
    g <- graph_from_incidence_matrix(cor.matrix, weighted=TRUE)
    estimated.celltypes <- as.vector(max_bipartite_match(g)$matching[components])
    # Assign
    colnames(W_RNA) <- estimated.celltypes
    rownames(H_RNA) <- estimated.celltypes
    rownames(H_Epi) <- estimated.celltypes
    list(W_RNA=W_RNA, H_RNA=H_RNA, H_Epi=H_Epi)
}

.assignCelltypeNames_List <- function(X_RNA, label, W_RNA, H_RNA, H_Epi){
    X_RNA2 <- do.call("rbind", X_RNA)
    W_RNA2 <- do.call("rbind", W_RNA)
    asn <- .assignCelltypeNames_Matrix(X_RNA2, label, W_RNA2, H_RNA, H_Epi)
    list(W_RNA=.Mat2List(X_RNA, asn$W_RNA),
        H_RNA=asn$H_RNA, H_Epi=asn$H_Epi)
}
