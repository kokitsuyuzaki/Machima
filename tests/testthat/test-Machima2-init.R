#
# Helper
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
# Test: init_W_RNA, init_H_RNA, init_H_Sym in single matrix mode
#
W0 <- matrix(runif(20*J), nrow=20, ncol=J)
H0 <- matrix(runif(J*30), nrow=J, ncol=30)
HS0 <- .makeSymMatrix(J)

out_init <- Machima2(X_RNA, X_Epi,
    init_W_RNA=W0, init_H_RNA=H0, init_H_Sym=HS0,
    J=J, num.iter=10)

expect_true(is.list(out_init))
expect_equal(dim(out_init$W_RNA), c(20, J))
expect_equal(dim(out_init$H_RNA), c(J, 30))
expect_equal(dim(out_init$H_Sym), c(J, J))
expect_true(isSymmetric(out_init$H_Sym))

#
# Test: all fix + all init => no updates, factors unchanged
#
out_fix <- Machima2(X_RNA, X_Epi,
    T=matrix(runif(15*20), 15, 20),
    init_W_RNA=W0, init_H_RNA=H0, init_H_Sym=HS0,
    fixW_RNA=TRUE, fixH_RNA=TRUE, fixT=TRUE, fixH_Sym=TRUE,
    J=J, num.iter=5)

expect_equal(out_fix$W_RNA, W0)
expect_equal(out_fix$H_RNA, H0)
expect_equal(out_fix$H_Sym, HS0)

#
# Test: two-stage workflow (Stage A: NMF, Stage B: fix RNA, fit H_Sym)
#
nmf_res <- nnTensor::NMF(X_RNA, J=J, num.iter=30, algorithm="Frobenius")
W_nmf <- nmf_res$U
H_nmf <- t(nmf_res$V)

out_2stage <- Machima2(X_RNA, X_Epi,
    init_W_RNA=W_nmf, init_H_RNA=H_nmf,
    fixW_RNA=TRUE, fixH_RNA=TRUE,
    J=J, num.iter=30)

# W_RNA should be exactly the NMF result (frozen)
expect_equal(out_2stage$W_RNA, W_nmf)
expect_equal(out_2stage$H_RNA, H_nmf)
# H_Sym should be fitted (not the random init)
expect_true(isSymmetric(out_2stage$H_Sym))

#
# Test: init in list mode
#
X_RNAs <- list(
    matrix(runif(20*30), nrow=20, ncol=30),
    matrix(runif(25*30), nrow=25, ncol=30))
X_Epis <- list(.makeSymMatrix(15), .makeSymMatrix(18))

W0_list <- list(
    matrix(runif(20*J), nrow=20, ncol=J),
    matrix(runif(25*J), nrow=25, ncol=J))
H0_list <- matrix(runif(J*30), nrow=J, ncol=30)
HS0_list <- .makeSymMatrix(J)

out_list <- Machima2(X_RNAs, X_Epis,
    init_W_RNA=W0_list, init_H_RNA=H0_list, init_H_Sym=HS0_list,
    fixW_RNA=TRUE, fixH_RNA=TRUE, fixH_Sym=TRUE,
    J=J, num.iter=5)

expect_equal(out_list$W_RNA, W0_list)
expect_equal(out_list$H_RNA, H0_list)
expect_equal(out_list$H_Sym, HS0_list)

#
# Test: NMF init params (n_restart, num_iter, algorithm)
#
out_nmf_params <- Machima2(X_RNA, X_Epi,
    init="RandomEpi",
    nmf_init_n_restart=2, nmf_init_num_iter=10,
    nmf_init_algorithm="Frobenius",
    J=J, num.iter=5)

expect_true(is.list(out_nmf_params))
expect_equal(dim(out_nmf_params$W_RNA), c(20, J))

#
# Test: Mat2List export
#
W_big <- matrix(runif(45*J), 45, J)
W_split <- Mat2List(X_RNAs, W_big)
expect_equal(length(W_split), 2)
expect_equal(nrow(W_split[[1]]), 20)
expect_equal(nrow(W_split[[2]]), 25)

#
# Test: fixH_Sym alone
#
out_fixH <- Machima2(X_RNA, X_Epi,
    init_H_Sym=HS0, fixH_Sym=TRUE,
    J=J, num.iter=10)
expect_equal(out_fixH$H_Sym, HS0)
