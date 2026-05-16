.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}
set.seed(42)
X_RNA <- matrix(runif(15*30), 15, 30)
X_Epi <- .makeSymMatrix(15)
J <- 3

# Default (use_shared_background=FALSE) reproduces v1.6.1
set.seed(1)
out_default <- Machima2(X_RNA, X_Epi, J=J, num.iter=5)
expect_null(out_default$w_0)
expect_null(out_default$delta)

# use_shared_background=TRUE: basic shape test
out_sh <- Machima2(X_RNA, X_Epi, J=J, num.iter=10,
    use_shared_background=TRUE)
expect_true(!is.null(out_sh$w_0))
expect_true(!is.null(out_sh$delta))
expect_equal(length(out_sh$w_0), 15)  # n_k vector
expect_equal(dim(out_sh$delta), c(15, J))
expect_equal(length(out_sh$h_vec), J + 1)
expect_true(all(out_sh$w_0 >= 0))
expect_true(all(out_sh$delta >= 0))
expect_true(all(out_sh$h_vec >= 0))

# Prediction functions
R_c1 <- out_sh$predict_celltype_R(1)
expect_equal(dim(R_c1), c(15, 15))
D_c1 <- out_sh$predict_celltype_delta(1)
expect_equal(dim(D_c1), c(15, 15))

# RecError decreases
errs <- out_sh$RecError[!is.na(out_sh$RecError)]
expect_true(errs[length(errs)] < errs[2])

# lambda_delta=Inf: delta not updated (stays at init)
set.seed(99)
out_inf1 <- Machima2(X_RNA, X_Epi, J=J, num.iter=1,
    use_shared_background=TRUE, lambda_delta=Inf)
set.seed(99)
out_inf5 <- Machima2(X_RNA, X_Epi, J=J, num.iter=5,
    use_shared_background=TRUE, lambda_delta=Inf)
# delta should be identical (not updated, same init)
expect_equal(out_inf1$delta, out_inf5$delta)

# List mode
X_RNAs <- list(matrix(runif(15*30), 15, 30), matrix(runif(18*30), 18, 30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))
out_list <- Machima2(X_RNAs, X_Epis, J=J, num.iter=5,
    use_shared_background=TRUE)
expect_equal(length(out_list$w_0), 2)
expect_equal(length(out_list$w_0[[1]]), 15)
expect_equal(length(out_list$w_0[[2]]), 18)
expect_equal(nrow(out_list$delta[[1]]), 15)
expect_equal(nrow(out_list$delta[[2]]), 18)

# Mutual exclusion with lambda_coupling
expect_error(
    Machima2(X_RNA, X_Epi, J=J,
        use_shared_background=TRUE, lambda_coupling=1),
    "mutually exclusive"
)
