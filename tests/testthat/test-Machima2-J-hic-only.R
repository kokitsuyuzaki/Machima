.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

set.seed(42)
X_RNA <- matrix(runif(15*30), 15, 30)
X_Epi <- .makeSymMatrix(15)
J <- 3

# J_hic_only=0 default
set.seed(1)
out0 <- Machima2(X_RNA, X_Epi, J=J, num.iter=10)
set.seed(1)
out0b <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, J_hic_only=0L)
expect_equal(out0$W_RNA, out0b$W_RNA)
expect_null(out0$W_hic)

# J_hic_only=2
out2 <- Machima2(X_RNA, X_Epi, J=J, num.iter=20, J_hic_only=2L)
expect_true(!is.null(out2$W_hic))
expect_equal(ncol(out2$W_hic), 2)
expect_equal(length(out2$h_hic), 2)
expect_true(all(out2$W_hic >= 0))
expect_true(all(out2$h_hic >= 0))
# H_Sym is still J x J
expect_equal(dim(out2$H_Sym), c(J, J))

# J_hic_only reduces RecError vs J_hic_only=0
set.seed(1)
out_no <- Machima2(X_RNA, X_Epi, J=J, num.iter=30, J_hic_only=0L)
set.seed(1)
out_yes <- Machima2(X_RNA, X_Epi, J=J, num.iter=30, J_hic_only=2L)
errs_no <- tail(out_no$RecError[!is.na(out_no$RecError)], 1)
errs_yes <- tail(out_yes$RecError[!is.na(out_yes$RecError)], 1)
expect_true(errs_yes <= errs_no)

# fixW_hic
W_hic_init <- matrix(runif(15*2, 0.1, 1), 15, 2)
out_fix <- Machima2(X_RNA, X_Epi, J=J, num.iter=10,
    J_hic_only=2L, W_hic_init=W_hic_init, fixW_hic=TRUE)
expect_equal(out_fix$W_hic, W_hic_init)

# List mode
X_RNAs <- list(matrix(runif(15*30), 15, 30), matrix(runif(18*30), 18, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_list <- Machima2(X_RNAs, X_Epis, J=J, num.iter=10, J_hic_only=2L)
expect_equal(length(out_list$W_hic), 2)
expect_equal(ncol(out_list$W_hic[[1]]), 2)
expect_equal(ncol(out_list$W_hic[[2]]), 2)
expect_equal(nrow(out_list$W_hic[[1]]), 15)
expect_equal(nrow(out_list$W_hic[[2]]), 18)
