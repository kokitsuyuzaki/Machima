######## 400 bins × 100 genes ########
set.seed(1234)
T = matrix(1E-10*runif(400*100), nrow=400, ncol=100)
# Block1
T[1:100, 1:30] <- 10*runif(100*30)
# Block2
T[201:300, 31:60] <- 10*runif(100*30)
# Block3
T[101:200, 61:90] <- 10*runif(100*30)

######## 100 genes × 300 cells ########
# Block1 (lamda100): [1:30, 1:90]
# Block2 (lamda100): [31:60, 91:180]
# Block3 (lamda100): [61:90, 181:270]
X_RNA = nnTensor::toyModel(model="siNMF_Easy")[[1]]
X_RNA[1:30, 1:90] <- rpois(30*90, lambda=140) # 140
X_RNA[31:60, 91:180] <- rpois(30*90, lambda=120) # 120

######## 100 genes × 20 samples ########
# Block1 (lamda100): [1:30, 61:120]
# Block2 (lamda100): [31:60, 121:180]
# Block3 (lamda100): [61:90, 1:60]
X_GAM = nnTensor::toyModel(model="siNMF_Easy")[[2]]
X_GAM[1:30, 61:120] <- rpois(30*60, lambda=140) # 140
X_GAM[31:60, 121:180] <- rpois(30*60, lambda=120) # 120
X_GAM = X_GAM[, c(1:5, 61:65, 121:125, 181:185)]

######## 400 bins × 20 samples ########
X_Epi <- T %*% X_GAM

#
# Single Matrix mode
#
out1 <- Machima(X_RNA, X_Epi, Beta=2, viz=TRUE)

#
# Single Matrix mode with T
#
out2 <- Machima(X_RNA, X_Epi, Beta=2, T=T, fixT=TRUE, viz=TRUE)

#
# Multiple Matrices mode
#
idx_RNA <- split(seq(nrow(X_RNA)), 1:3)
idx_Epi <- split(seq(nrow(X_Epi)), 1:3)
X_RNAs <- lapply(idx_RNA, function(x){X_RNA[x, ]})
X_Epis <- lapply(idx_Epi, function(x){X_Epi[x, ]})

out3 <- Machima(X_RNAs, X_Epis, Beta=2, viz=TRUE)

#
# Multiple Matrices mode with T
#
Ts <- lapply(seq_along(idx_RNA), function(x){T[idx_Epi[[x]], idx_RNA[[x]]]})
out4 <- Machima(X_RNAs, X_Epis, T=Ts, fixT=TRUE, Beta=2, viz=TRUE)
