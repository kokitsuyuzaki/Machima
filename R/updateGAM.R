.updateGAM <- function(X_RNA, X_Epi, W_RNA, H_Epi, T, horizontal){
    if(horizontal){
        .updateGAM_HZL(X_RNA, X_Epi, W_RNA, H_Epi, T)
    }else{
        .updateGAM_Normal(X_RNA, X_Epi, W_RNA, H_Epi, T)
    }
}

.updateGAM_HZL <- function(X_RNA, X_Epi, W_RNA, H_Epi, T){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        X_GAM <- t(T) %*% X_Epi
    }else{
        X_GAM <- lapply(seq_along(X_Epi), function(x){
            t(T[[x]]) %*% X_Epi[[x]]
        })
    }
    X_GAM
}

.updateGAM_Normal <- function(X_RNA, X_Epi, W_RNA, H_Epi, T){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        X_GAM <- W_RNA %*% H_Epi
    }else{
        W_RNA2 <- do.call("rbind", W_RNA)
        X_GAM <- W_RNA2 %*% H_Epi
    }
    X_GAM
}
