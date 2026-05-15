# Helper to compute full Hi-C reconstruction including hic-only contribution
# R_full = G · H_Sym · Gᵀ + W_hic · diag(h_hic) · W_hicᵀ

.addHic <- function(S_hat, W_hic_k, h_hic){
    if(length(h_hic) > 0L && !is.null(W_hic_k)){
        S_hat + W_hic_k %*% (h_hic * t(W_hic_k))
    }else{
        S_hat
    }
}

.reconstructEpi_single <- function(W_RNA_k, T_k, W_hic_k, H_Sym, h_hic, U_k=NULL){
    W_E <- if(!is.null(U_k)) W_RNA_k + U_k else W_RNA_k
    G <- T_k %*% W_E
    R <- G %*% H_Sym %*% t(G)
    .addHic(R, W_hic_k, h_hic)
}
