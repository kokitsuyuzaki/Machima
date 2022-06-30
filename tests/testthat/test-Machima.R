#
# Single Matrix mode
#
X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
X_Epi <- matrix(runif(15*25), nrow=15, ncol=25)

out1 <- Machima(X_RNA, X_Epi, T=NULL)

expect_true(is.list(out1))
expect_equal(dim(out1$W_RNA), c(20, 3))
expect_equal(dim(out1$H_RNA), c(3, 30))
expect_equal(dim(out1$H_Epi), c(3, 25))
expect_equal(dim(out1$X_GAM), c(20, 25))
expect_equal(dim(out1$T), c(15, 20))

#
# Single Matrix mode with T
#
T <- matrix(runif(15*20), nrow=15, ncol=20)

out2 <- Machima(X_RNA, X_Epi, T=T)

expect_true(is.list(out2))
expect_equal(dim(out2$W_RNA), c(20, 3))
expect_equal(dim(out2$H_RNA), c(3, 30))
expect_equal(dim(out2$H_Epi), c(3, 25))
expect_equal(dim(out2$X_GAM), c(20, 25))
expect_equal(dim(out2$T), c(15, 20))


#
# Multiple Matrices mode
#
X_RNAs <- list(
    X_RNA_1 = matrix(runif(20*30), nrow=20, ncol=30),
    X_RNA_2 = matrix(runif(30*30), nrow=30, ncol=30),
    X_RNA_3 = matrix(runif(25*30), nrow=25, ncol=30),
    X_RNA_4 = matrix(runif(35*30), nrow=35, ncol=30),
    X_RNA_5 = matrix(runif(40*30), nrow=40, ncol=30))

X_Epis <- list(
    X_Epi_1 = matrix(runif(25*25), nrow=25, ncol=25),
    X_Epi_2 = matrix(runif(35*25), nrow=35, ncol=25),
    X_Epi_3 = matrix(runif(20*25), nrow=20, ncol=25),
    X_Epi_4 = matrix(runif(30*25), nrow=30, ncol=25),
    X_Epi_5 = matrix(runif(40*25), nrow=40, ncol=25))

out3 <- Machima(X_RNAs, X_Epis, T=NULL)

expect_true(is.list(out3))
expect_true(is.list(out3$W_RNA))
expect_equal(dim(out3$W_RNA[[1]]), c(20, 3))
expect_equal(dim(out3$W_RNA[[2]]), c(30, 3))
expect_equal(dim(out3$W_RNA[[3]]), c(25, 3))
expect_equal(dim(out3$W_RNA[[4]]), c(35, 3))
expect_equal(dim(out3$W_RNA[[5]]), c(40, 3))
expect_equal(dim(out3$H_RNA), c(3, 30))
expect_equal(dim(out3$H_Epi), c(3, 25))
expect_equal(dim(out3$X_GAM), c(150, 25))
expect_true(is.list(out3$T))
expect_equal(dim(out3$T[[1]]), c(25, 20))
expect_equal(dim(out3$T[[2]]), c(35, 30))
expect_equal(dim(out3$T[[3]]), c(20, 25))
expect_equal(dim(out3$T[[4]]), c(30, 35))
expect_equal(dim(out3$T[[5]]), c(40, 40))

#
# Multiple Matrices mode
#
Ts <- list(
	T_1 = matrix(runif(25*20), nrow=25, ncol=20),
	T_2 = matrix(runif(35*30), nrow=35, ncol=30),
	T_3 = matrix(runif(20*25), nrow=20, ncol=25),
	T_4 = matrix(runif(30*35), nrow=30, ncol=35),
	T_5 = matrix(runif(40*40), nrow=40, ncol=40))

out4 <- Machima(X_RNAs, X_Epis, T=Ts)

expect_true(is.list(out4))
expect_true(is.list(out4$W_RNA))
expect_equal(dim(out4$W_RNA[[1]]), c(20, 3))
expect_equal(dim(out4$W_RNA[[2]]), c(30, 3))
expect_equal(dim(out4$W_RNA[[3]]), c(25, 3))
expect_equal(dim(out4$W_RNA[[4]]), c(35, 3))
expect_equal(dim(out4$W_RNA[[5]]), c(40, 3))
expect_equal(dim(out4$H_RNA), c(3, 30))
expect_equal(dim(out4$H_Epi), c(3, 25))
expect_equal(dim(out4$X_GAM), c(150, 25))
expect_true(is.list(out4$T))
expect_equal(dim(out4$T[[1]]), c(25, 20))
expect_equal(dim(out4$T[[2]]), c(35, 30))
expect_equal(dim(out4$T[[3]]), c(20, 25))
expect_equal(dim(out4$T[[4]]), c(30, 35))
expect_equal(dim(out4$T[[5]]), c(40, 40))
