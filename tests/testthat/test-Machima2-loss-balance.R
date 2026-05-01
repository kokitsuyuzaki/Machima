#
# Test lambda_balance
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
# Default equivalence: lambda_balance=0.5 == no arg
#
set.seed(1)
out_default <- Machima2(X_RNA, X_Epi, J=J, num.iter=10)
set.seed(1)
out_half <- Machima2(X_RNA, X_Epi, J=J, num.iter=10, lambda_balance=0.5)
expect_equal(out_default$W_RNA, out_half$W_RNA)
expect_equal(out_default$H_RNA, out_half$H_RNA)
expect_equal(out_default$H_Sym, out_half$H_Sym)

#
# Validation
#
expect_error(Machima2(X_RNA, X_Epi, J=J, lambda_balance=-0.1))
expect_error(Machima2(X_RNA, X_Epi, J=J, lambda_balance=1.5))

#
# lambda_balance works in Machima too
#
X_Epi2 <- matrix(runif(15*25), 15, 25)
set.seed(1)
out_m_default <- Machima(X_RNA, X_Epi2, num.iter=5)
set.seed(1)
out_m_half <- Machima(X_RNA, X_Epi2, num.iter=5, lambda_balance=0.5)
expect_equal(out_m_default$W_RNA, out_m_half$W_RNA)
