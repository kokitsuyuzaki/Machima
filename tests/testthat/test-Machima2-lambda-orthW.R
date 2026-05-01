#
# Test lambda_orthW (continuous orthogonality strength)
#
.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

set.seed(42)
X_RNA <- matrix(runif(15*30), nrow=15, ncol=30)
X_Epi <- .makeSymMatrix(15)
J <- 3

#
# lambda_orthW=0 reproduces default (no orthogonality)
#
set.seed(1)
out_default <- Machima2(X_RNA, X_Epi, J=J, num.iter=10)
set.seed(1)
out_zero <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, lambda_orthW=0)
expect_equal(out_default$W_RNA, out_zero$W_RNA)

#
# lambda_orthW=1 reproduces orthW_RNA=TRUE
#
set.seed(1)
suppressWarnings(out_orth <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, orthW_RNA=TRUE))
set.seed(1)
out_one <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, lambda_orthW=1)
expect_equal(out_orth$W_RNA, out_one$W_RNA)

#
# orthW_RNA=TRUE triggers deprecation warning
#
expect_warning(
    Machima2(X_RNA, X_Epi, J=J, num.iter=1, orthW_RNA=TRUE),
    "deprecated"
)

#
# Intermediate values produce valid output
#
out_mid <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, lambda_orthW=0.3)
expect_true(all(out_mid$W_RNA >= 0))
expect_false(any(is.nan(out_mid$W_RNA)))

#
# Argument validation
#
expect_error(Machima2(X_RNA, X_Epi, J=J, lambda_orthW=-0.1))
expect_error(Machima2(X_RNA, X_Epi, J=J, lambda_orthW=1.5))

#
# Works in Machima (asymmetric) too
#
X_Epi2 <- matrix(runif(15*25), 15, 25)
expect_warning(
    Machima(X_RNA, X_Epi2, num.iter=1, orthW_RNA=TRUE),
    "deprecated"
)
out_m <- Machima(X_RNA, X_Epi2, num.iter=5, lambda_orthW=0.5)
expect_true(all(out_m$W_RNA >= 0))
