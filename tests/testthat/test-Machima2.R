#
# Helper to create symmetric non-negative matrix
#
.makeSymMatrix <- function(n){
    A <- matrix(runif(n*n), nrow=n, ncol=n)
    (A + t(A)) / 2
}

#
# Single Matrix mode
#
X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
X_Epi <- .makeSymMatrix(15)

out1 <- Machima2(X_RNA, X_Epi)

expect_true(is.list(out1))
expect_equal(dim(out1$W_RNA), c(20, 3))
expect_equal(dim(out1$H_RNA), c(3, 30))
expect_equal(dim(out1$H_Sym), c(3, 3))
expect_true(isSymmetric(out1$H_Sym))
expect_equal(dim(out1$T), c(15, 20))

#
# Single Matrix mode with T
#
T <- matrix(runif(15*20), nrow=15, ncol=20)

out2 <- Machima2(X_RNA, X_Epi, T=T)

expect_true(is.list(out2))
expect_equal(dim(out2$W_RNA), c(20, 3))
expect_equal(dim(out2$H_RNA), c(3, 30))
expect_equal(dim(out2$H_Sym), c(3, 3))
expect_true(isSymmetric(out2$H_Sym))
expect_equal(dim(out2$T), c(15, 20))

#
# Single Matrix mode with horizontal
#
out2h <- Machima2(X_RNA, X_Epi, T=T, horizontal=TRUE)

expect_true(is.list(out2h))
expect_equal(dim(out2h$W_RNA), c(20, 3))
expect_equal(dim(out2h$H_RNA), c(3, 30))
expect_equal(dim(out2h$H_Sym), c(3, 3))
expect_true(isSymmetric(out2h$H_Sym))

#
# Multiple Matrices mode
#
X_RNAs <- list(
    X_RNA_1 = matrix(runif(20*30), nrow=20, ncol=30),
    X_RNA_2 = matrix(runif(30*30), nrow=30, ncol=30),
    X_RNA_3 = matrix(runif(25*30), nrow=25, ncol=30))

X_Epis <- list(
    X_Epi_1 = .makeSymMatrix(25),
    X_Epi_2 = .makeSymMatrix(35),
    X_Epi_3 = .makeSymMatrix(20))

out3 <- Machima2(X_RNAs, X_Epis, T=NULL)

expect_true(is.list(out3))
expect_true(is.list(out3$W_RNA))
expect_equal(dim(out3$W_RNA[[1]]), c(20, 3))
expect_equal(dim(out3$W_RNA[[2]]), c(30, 3))
expect_equal(dim(out3$W_RNA[[3]]), c(25, 3))
expect_equal(dim(out3$H_RNA), c(3, 30))
expect_equal(dim(out3$H_Sym), c(3, 3))
expect_true(isSymmetric(out3$H_Sym))
expect_true(is.list(out3$T))
expect_equal(dim(out3$T[[1]]), c(25, 20))
expect_equal(dim(out3$T[[2]]), c(35, 30))
expect_equal(dim(out3$T[[3]]), c(20, 25))

#
# Multiple Matrices mode with Ts
#
Ts <- list(
    T_1 = matrix(runif(25*20), nrow=25, ncol=20),
    T_2 = matrix(runif(35*30), nrow=35, ncol=30),
    T_3 = matrix(runif(20*25), nrow=20, ncol=25))

out4 <- Machima2(X_RNAs, X_Epis, T=Ts)

expect_true(is.list(out4))
expect_true(is.list(out4$W_RNA))
expect_equal(dim(out4$W_RNA[[1]]), c(20, 3))
expect_equal(dim(out4$W_RNA[[2]]), c(30, 3))
expect_equal(dim(out4$W_RNA[[3]]), c(25, 3))
expect_equal(dim(out4$H_RNA), c(3, 30))
expect_equal(dim(out4$H_Sym), c(3, 3))
expect_true(isSymmetric(out4$H_Sym))
expect_true(is.list(out4$T))
expect_equal(dim(out4$T[[1]]), c(25, 20))
expect_equal(dim(out4$T[[2]]), c(35, 30))
expect_equal(dim(out4$T[[3]]), c(20, 25))

#
# Reconstruction error decreases
#
out5 <- Machima2(X_RNA, X_Epi, num.iter=50)
errors <- out5$RecError[!is.na(out5$RecError)]
expect_true(length(errors) > 1)
