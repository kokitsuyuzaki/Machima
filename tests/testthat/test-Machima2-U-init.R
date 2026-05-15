.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}
set.seed(42)
X_RNA <- matrix(runif(15*30), 15, 30)
X_Epi <- .makeSymMatrix(15)
J <- 3

# lambda_coupling=Inf: bit-identical to v1.5.0 (no U consumed)
set.seed(1)
out_inf <- Machima2(X_RNA, X_Epi, J=J, num.iter=5)
expect_null(out_inf$U)

# lambda_coupling=1, init_U=NULL: U should be non-zero after iteration
out_soft <- Machima2(X_RNA, X_Epi, J=J, num.iter=5, lambda_coupling=1)
expect_true(all(out_soft$U > 0))

# Reproducibility: same seed → same U
set.seed(10)
out_a <- Machima2(X_RNA, X_Epi, J=J, num.iter=3, lambda_coupling=1)
set.seed(10)
out_b <- Machima2(X_RNA, X_Epi, J=J, num.iter=3, lambda_coupling=1)
expect_equal(out_a$U, out_b$U)

# Different seeds → different U
set.seed(10)
out_c <- Machima2(X_RNA, X_Epi, J=J, num.iter=3, lambda_coupling=1)
set.seed(20)
out_d <- Machima2(X_RNA, X_Epi, J=J, num.iter=3, lambda_coupling=1)
expect_false(identical(out_c$U, out_d$U))

# Explicit init_U overrides the auto-default
U_user <- matrix(runif(15*J, 0.5, 1), 15, J)
out_explicit <- Machima2(X_RNA, X_Epi, J=J, num.iter=3,
    lambda_coupling=1, init_U=U_user, fixU=TRUE)
expect_equal(out_explicit$U, U_user)
