X_RNA <- matrix(runif(20*30), nrow=20, ncol=30)
X_Epi <- matrix(runif(15*25), nrow=15, ncol=25)

out <- Machima(X_RNA, X_Epi, T=NULL, verbose=TRUE)

expect_true(is.list(out))