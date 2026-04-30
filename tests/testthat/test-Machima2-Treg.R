#
# Test T_regularization options
#
.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

set.seed(123)
X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
X_Epi <- .makeSymMatrix(15)
J <- 3

#
# frobenius_unit: H_Sym should have meaningful magnitude
#
out_fu <- Machima2(X_RNA, X_Epi,
    T_regularization="frobenius_unit",
    J=J, num.iter=30)

expect_true(is.list(out_fu))
expect_equal(dim(out_fu$H_Sym), c(J, J))
expect_true(isSymmetric(out_fu$H_Sym))
# T should be roughly unit Frobenius norm
if(is.matrix(out_fu$T)){
    expect_true(abs(norm(out_fu$T, "F") - 1) < 0.1)
}
# RecError should still decrease
errs <- out_fu$RecError[!is.na(out_fu$RecError)]
expect_true(errs[length(errs)] < errs[2])

#
# l2: lambda_T penalty
#
out_l2 <- Machima2(X_RNA, X_Epi,
    T_regularization="l2", lambda_T=0.1,
    J=J, num.iter=30)

expect_true(is.list(out_l2))
expect_true(isSymmetric(out_l2$H_Sym))

#
# none: backward compatible (default)
#
out_none <- Machima2(X_RNA, X_Epi,
    T_regularization="none",
    J=J, num.iter=10)
expect_true(is.list(out_none))

#
# frobenius_unit with fixT=TRUE should be ignored (T unchanged)
#
T_fixed <- matrix(runif(15*20), 15, 20)
out_fixT <- Machima2(X_RNA, X_Epi,
    T=T_fixed, fixT=TRUE,
    T_regularization="frobenius_unit",
    J=J, num.iter=5)
expect_equal(out_fixT$T, T_fixed)

#
# List mode with frobenius_unit
#
X_RNAs <- list(matrix(runif(20*30), 20, 30), matrix(runif(25*30), 25, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_list_fu <- Machima2(X_RNAs, X_Epis,
    T_regularization="frobenius_unit",
    J=J, num.iter=20)

expect_true(is.list(out_list_fu))
expect_true(isSymmetric(out_list_fu$H_Sym))
# Each T should be roughly unit norm
for(i in seq_along(out_list_fu$T)){
    expect_true(abs(norm(out_list_fu$T[[i]], "F") - 1) < 0.1)
}
