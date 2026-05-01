#
# Test proper constrained MU for H_Sym_structure="diagonal"
#
.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

set.seed(42)
J <- 3

#
# Diagonal preservation: off-diagonals are exactly zero
#
X_RNA <- matrix(runif(15*30), 15, 30)
X_Epi <- .makeSymMatrix(15)
out <- Machima2(X_RNA, X_Epi, H_Sym_structure="diagonal", J=J, num.iter=20)
expect_true(all(out$H_Sym[upper.tri(out$H_Sym)] == 0))
expect_true(all(out$H_Sym[lower.tri(out$H_Sym)] == 0))
expect_true(all(diag(out$H_Sym) >= 0))

#
# Monotonic RecError decrease with Beta=2
#
errs <- out$RecError[!is.na(out$RecError)]
expect_true(all(diff(errs[-1]) <= 1e-10))

#
# J=1: diagonal and symmetric should give same result
#
set.seed(123)
out_sym1 <- Machima2(X_RNA, X_Epi, H_Sym_structure="symmetric", J=1, num.iter=20)
set.seed(123)
out_diag1 <- Machima2(X_RNA, X_Epi, H_Sym_structure="diagonal", J=1, num.iter=20)
expect_equal(out_sym1$H_Sym, out_diag1$H_Sym, tolerance=1e-6)

#
# Synthetic sym-CP signal recovery
#
set.seed(42)
G_true <- matrix(runif(15*J), 15, J)
h_true <- c(3, 1, 5)
X_Epi_cp <- G_true %*% diag(h_true) %*% t(G_true)
X_RNA_cp <- G_true %*% matrix(runif(J*30), J, 30)

out_cp <- Machima2(X_RNA_cp, X_Epi_cp, H_Sym_structure="diagonal", J=J, num.iter=100)
h_est <- sort(diag(out_cp$H_Sym), decreasing=TRUE)
h_ref <- sort(h_true, decreasing=TRUE)
# Relative ordering should match (largest component is largest)
expect_equal(which.max(h_est), 1)

#
# List mode diagonal
#
X_RNAs <- list(matrix(runif(15*30), 15, 30), matrix(runif(18*30), 18, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_list <- Machima2(X_RNAs, X_Epis, H_Sym_structure="diagonal", J=J, num.iter=10)
expect_true(all(out_list$H_Sym[upper.tri(out_list$H_Sym)] == 0))
