# Update GAM matrix
.updateGAM <- function(X_RNA, X_Epi, W_RNA, H_Epi){
    if(is.matrix(X_RNA) && is.matrix(X_Epi)){
        X_GAM <- W_RNA %*% H_Epi
    }else{
        W_RNA2 <- do.call("rbind", W_RNA)
        X_GAM <- W_RNA2 %*% H_Epi
    }
    X_GAM
}
