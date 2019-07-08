rm(list = ls())


load("data/binaClinical.RData")
load("data/binaMutation.RData")

# SVD
(s <- svd(t(binaMutation)) ) #genes x patients

# Get diagonal
D <- diag(s$d)

# Getting the top 50% of the diagonal 
Pdim = ceiling(dim(s$u)[2]*0.5)

# Matrix recunstriction via top diag multiplication
X_smooth = s$u[,1:Pdim] %*% D[1:Pdim,1:Pdim] %*% t(s$v)[1:Pdim,] #  X = U D V'

AbsoluteErrorAproxSum = (sum(binaMutation)-sum(X_smooth))
RelativeErrorAproxSum = (sum(binaMutation)-sum(X_smooth))/sum(binaMutation)
cat("\n\nRelativeErrorAproxSum  :", RelativeErrorAproxSum)

matrixbM = t(as.matrix(binaMutation))
row.names(matrixbM) <- NULL
names(matrixbM) <- NULL
AbsoluteErrorAproxNorm = max(apply(abs(matrixbM - X_smooth ),1,sum))
RelativeErrorAproxNorm = max(apply(abs(matrixbM - X_smooth ),1,sum)) / max(apply(abs(matrixbM),1,sum))
cat("\n\nRelativeErrorAproxNorm  :", RelativeErrorAproxNorm)

mutationSmooth <- as.data.frame(t(X_smooth))
row.names(mutationSmooth) <- row.names(binaMutation)
names(mutationSmooth) <- names(binaMutation)
save(mutationSmooth, file="data/mutationSmooth.RData")

#write csv file for python script
write.table(X_smooth, file = "temp/10-binaMsmooth.csv", sep =';', dec = '.', row.names = FALSE )
