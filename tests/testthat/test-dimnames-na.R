#
# Test: dimnames never contain NA after label assignment
#
set.seed(42)

# Single matrix: more components than cell types (J > num_celltypes)
X_RNA <- matrix(runif(30*50), nrow=30, ncol=50)
S <- matrix(runif(20*20), 20, 20)
X_Epi <- (S + t(S)) / 2
label <- rep(c("TypeA", "TypeB"), each=25)

# J=5 > 2 cell types: some components won't match any cell type
out <- Machima2(X_RNA, X_Epi, label=label, J=5, num.iter=10)

expect_false(any(is.na(rownames(out$H_RNA))))
expect_false(any(is.na(rownames(out$H_Sym))))
expect_false(any(is.na(colnames(out$H_Sym))))
expect_false(any(is.na(colnames(out$W_RNA))))

# Same test for original Machima
X_Epi2 <- matrix(runif(20*15), 20, 15)
out2 <- Machima(X_RNA, X_Epi2, label=label, J=5, num.iter=10)
expect_false(any(is.na(rownames(out2$H_RNA))))
expect_false(any(is.na(rownames(out2$H_Epi))))
expect_false(any(is.na(colnames(out2$W_RNA))))

# List mode
X_RNAs <- list(matrix(runif(20*50), 20, 50), matrix(runif(25*50), 25, 50))
X_Epis <- list((function(){A<-matrix(runif(15*15),15,15);(A+t(A))/2})(),
               (function(){A<-matrix(runif(18*18),18,18);(A+t(A))/2})())
out3 <- Machima2(X_RNAs, X_Epis, label=label, J=5, num.iter=10)
expect_false(any(is.na(rownames(out3$H_RNA))))
expect_false(any(is.na(rownames(out3$H_Sym))))
expect_false(any(is.na(colnames(out3$H_Sym))))
