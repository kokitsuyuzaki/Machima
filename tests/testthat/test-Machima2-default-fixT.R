#
# Test fixT=TRUE default and auto-identity T
#
.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

set.seed(42)
J <- 3

#
# New default: fixT=TRUE, auto-identity when l==n
#
X_RNA_sq <- matrix(runif(15*30), nrow=15, ncol=30)
X_Epi_sq <- .makeSymMatrix(15)

out_default <- Machima2(X_RNA_sq, X_Epi_sq, J=J, num.iter=10)
# T should be identity (15x15)
expect_equal(out_default$T, diag(15))

#
# Auto-identity in list mode (l==n per chrom)
#
X_RNAs_sq <- list(matrix(runif(15*30), 15, 30), matrix(runif(18*30), 18, 30))
X_Epis_sq <- list(.makeSymMatrix(15), .makeSymMatrix(18))

out_list <- Machima2(X_RNAs_sq, X_Epis_sq, J=J, num.iter=5)
expect_equal(out_list$T[[1]], diag(15))
expect_equal(out_list$T[[2]], diag(18))

#
# l != n: T is auto-constructed as random (not identity), still fixed
#
X_RNA_rect <- matrix(runif(20*30), nrow=20, ncol=30)
X_Epi_rect <- .makeSymMatrix(15)

out_rect <- Machima2(X_RNA_rect, X_Epi_rect, J=J, num.iter=10)
expect_equal(dim(out_rect$T), c(15, 20))
# Should not be identity (different dims)

#
# Explicit fixT=FALSE: learned T, deprecation warning
#
expect_warning(
    Machima2(X_RNA_rect, X_Epi_rect, fixT=FALSE, J=J, num.iter=5),
    "fixT = FALSE"
)

#
# No warning when T_regularization is set
#
out_lr <- Machima2(X_RNA_rect, X_Epi_rect,
    fixT=FALSE, T_regularization="low_rank", T_rank=4,
    J=J, num.iter=5)
expect_true(is.list(out_lr))

#
# No warning when user supplies T explicitly
#
T_user <- matrix(runif(15*20), 15, 20)
out_userT <- Machima2(X_RNA_rect, X_Epi_rect,
    T=T_user, fixT=FALSE, J=J, num.iter=5)
expect_true(is.list(out_userT))
