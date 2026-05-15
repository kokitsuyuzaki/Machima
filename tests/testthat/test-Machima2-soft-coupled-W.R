.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}
set.seed(42)
X_RNA <- matrix(runif(15*30), 15, 30)
X_Epi <- .makeSymMatrix(15)
J <- 3

# Default (Inf) reproduces v1.5.0
set.seed(1)
out_inf <- Machima2(X_RNA, X_Epi, J=J, num.iter=10)
set.seed(1)
out_inf2 <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, lambda_coupling=Inf)
expect_equal(out_inf$W_RNA, out_inf2$W_RNA)
expect_null(out_inf$U)

# fixU=TRUE, explicit init_U=zero → same as hard share
set.seed(1)
out_fix0 <- Machima2(X_RNA, X_Epi, J=J, num.iter=10,
    lambda_coupling=1.0, fixU=TRUE,
    init_U=matrix(0, 15, J))
expect_equal(out_fix0$W_RNA, out_inf$W_RNA)

# lambda_coupling < Inf, fixU=FALSE → U learned
out_soft <- Machima2(X_RNA, X_Epi, J=J, num.iter=20,
    lambda_coupling=1.0)
expect_true(!is.null(out_soft$U))
expect_equal(dim(out_soft$U), c(15, J))
expect_true(all(out_soft$U >= 0))

# lambda_coupling=0 works (no NaN)
out_zero <- Machima2(X_RNA, X_Epi, J=J, num.iter=10,
    lambda_coupling=0)
expect_false(any(is.nan(out_zero$U)))

# List mode
X_RNAs <- list(matrix(runif(15*30), 15, 30), matrix(runif(18*30), 18, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_list <- Machima2(X_RNAs, X_Epis, J=J, num.iter=10,
    lambda_coupling=1.0)
expect_equal(length(out_list$U), 2)
expect_equal(nrow(out_list$U[[1]]), 15)
expect_equal(nrow(out_list$U[[2]]), 18)
