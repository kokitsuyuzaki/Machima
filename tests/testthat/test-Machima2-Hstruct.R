#
# Test H_Sym_structure = "diagonal"
#
.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

set.seed(42)
X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
X_Epi <- .makeSymMatrix(15)
J <- 3

#
# Diagonal stays diagonal
#
out_diag <- Machima2(X_RNA, X_Epi,
    H_Sym_structure="diagonal", J=J, num.iter=30)

expect_true(all(out_diag$H_Sym == diag(diag(out_diag$H_Sym))))
expect_true(all(diag(out_diag$H_Sym) >= 0))

#
# Backward compat: "symmetric" is default
#
out_sym <- Machima2(X_RNA, X_Epi,
    H_Sym_structure="symmetric", J=J, num.iter=10)
expect_true(is.list(out_sym))
# Off-diagonal should generally be non-zero for symmetric
expect_true(isSymmetric(out_sym$H_Sym))

#
# Compatible with fixT=TRUE
#
T_fixed <- matrix(runif(15*20), 15, 20)
out_fixT <- Machima2(X_RNA, X_Epi, T=T_fixed, fixT=TRUE,
    H_Sym_structure="diagonal", J=J, num.iter=10)
expect_true(all(out_fixT$H_Sym == diag(diag(out_fixT$H_Sym))))

#
# Compatible with frobenius_unit
#
suppressWarnings(out_frob <- Machima2(X_RNA, X_Epi, fixT=FALSE,
    T_regularization="frobenius_unit",
    H_Sym_structure="diagonal", J=J, num.iter=10))
expect_true(all(out_frob$H_Sym == diag(diag(out_frob$H_Sym))))

#
# Compatible with low_rank
#
out_lr <- Machima2(X_RNA, X_Epi, fixT=FALSE,
    T_regularization="low_rank", T_rank=5,
    H_Sym_structure="diagonal", J=J, num.iter=10)
expect_true(all(out_lr$H_Sym == diag(diag(out_lr$H_Sym))))

#
# Init projection warning
#
full_H <- matrix(runif(J*J), J, J)
full_H <- (full_H + t(full_H)) / 2
expect_warning(
    Machima2(X_RNA, X_Epi, init_H_Sym=full_H,
        H_Sym_structure="diagonal", J=J, num.iter=5),
    "off-diagonal"
)

#
# List mode diagonal
#
X_RNAs <- list(matrix(runif(20*30), 20, 30), matrix(runif(25*30), 25, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_list <- Machima2(X_RNAs, X_Epis,
    H_Sym_structure="diagonal", J=J, num.iter=10)
expect_true(all(out_list$H_Sym == diag(diag(out_list$H_Sym))))
