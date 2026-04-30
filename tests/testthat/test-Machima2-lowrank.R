#
# Test T_regularization = "low_rank"
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
# Basic low_rank run
#
out_lr <- Machima2(X_RNA, X_Epi,
    T_regularization="low_rank", T_rank=5,
    J=J, num.iter=20)

expect_true(is.list(out_lr))
expect_equal(dim(out_lr$H_Sym), c(J, J))
expect_true(isSymmetric(out_lr$H_Sym))
# T_factors returned
expect_true(!is.null(out_lr$T_factors))
expect_equal(ncol(out_lr$T_factors$U), 5)
expect_equal(ncol(out_lr$T_factors$V), 5)
# T rank bounded
expect_true(qr(out_lr$T)$rank <= 5)
# RecError should not have NaN
errs <- out_lr$RecError[!is.na(out_lr$RecError)]
expect_false(any(is.nan(errs)))

#
# List mode low_rank
#
X_RNAs <- list(matrix(runif(20*30), 20, 30), matrix(runif(25*30), 25, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_lr_list <- Machima2(X_RNAs, X_Epis,
    T_regularization="low_rank", T_rank=4,
    J=J, num.iter=10)

expect_true(is.list(out_lr_list$T_factors))
expect_equal(length(out_lr_list$T_factors$U), 2)
expect_equal(ncol(out_lr_list$T_factors$U[[1]]), 4)
for(k in 1:2){
    expect_true(qr(out_lr_list$T[[k]])$rank <= 4)
}

#
# fixT + low_rank should error
#
expect_error(
    Machima2(X_RNA, X_Epi, fixT=TRUE,
        T_regularization="low_rank", T_rank=3, J=J),
    "incompatible")

#
# T_rank=NULL with low_rank should error
#
expect_error(
    Machima2(X_RNA, X_Epi,
        T_regularization="low_rank", J=J),
    "T_rank")
