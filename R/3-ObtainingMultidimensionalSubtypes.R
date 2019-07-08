rm(list = ls())
k<-3

# LoadFile --------
load(paste("dataExample/Filtered/bcCPk",k,".Rd",sep=''), envir = globalenv())

#' Projection of the factor matrices
# projection -------
A <- bcCP$A  #patients x k
B <- bcCP$B  #genes x k
C <- bcCP$C  #clinical vars x k


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

modeF <- A
plot(modeF[, 1], type = "n", ylim = range(modeF),
     xlab = "Patient Factor", ylab = "Component score")
points(modeF[, 1], pch = 16, col="red")
points(modeF[, 2], pch = 17, col="blue")
points(modeF[, 3], pch = 18)
legend("bottomright", legend = c("C1", "C2", "C3"), pch = c(16, 17,18), 
       col = c("red","blue", "black"), inset=c(-0.2,0),
       bty = "o", title = "component")

modeF <- B
plot(modeF[, 1], type = "n", ylim = range(modeF),
     xlab = "Gene Factor", ylab = "Component score")
points(modeF[, 1], pch = 16, col="red")
points(modeF[, 2], pch = 17, col="blue")
points(modeF[, 3], pch = 18)
legend("bottomright", legend = c("C1", "C2", "C3"), pch = c(16, 17,18), 
       col = c("red","blue", "black"), inset=c(-0.2,0),
       bty = "o", title = "component")

modeF <- C
pdf( file = paste0("temp/C",k,".pdf"),  onefile = TRUE, width = 9, height = 7)
plot(modeF[, 1], type = "n", ylim = range(modeF),
     xlab = "Clinical variable Factor", ylab = "Component score")
points(modeF[, 1], pch = 16, col="red")
points(modeF[, 2], pch = 17, col="blue")
points(modeF[, 3], pch = 18)
legend("bottomright", legend = c("C1", "C2", "C3"), pch = c(16, 17,18), 
       col = c("red","blue", "black"), inset=c(-0.2,0),
       bty = "o", title = "component")
dev.off()




pdf( file = paste0("temp/allComponents",k,".pdf"),  onefile = TRUE, width = 9, height = 7)

#u = A
#v = B
#w = C
data_frame(component = c(rep("u[1]", 482), rep("u[2]", 482), rep("u[3]", 482),
                         rep("v[1]", 499), rep("v[2]", 499), rep("v[3]", 499),
                         rep("w[1]", 32), rep("w[2]", 32), rep("w[3]", 32)),
           value = c(A[ , 1], B[ , 1], C[ , 1],
                     A[ , 2], B[ , 2], C[ , 2],
                     A[ , 3], B[ , 3], C[ , 3]),
           index = c(rep(1:482, 3),
                     rep(1:499, 3), 
                     rep(1:32, 3))) %>%
  ggplot(aes(index, value)) + geom_line() +
  facet_wrap(~component, scales = "free", nrow = 3,
             labeller = labeller(component = label_parsed)) +
  theme(axis.title = element_blank())
dev.off()



#Plot all components ---------

rm(list = ls())
load("dataExample/Filtered/MatrixOrderPxGxC.RData")

X <- array(NA, dim = c(32, 482, 499))
indexG=1
for(i in (1:32) ) {
  X[i, , ]    <- array(unlist(MatrixOrderPxGxC[ , indexG:(indexG+498)])) 
  indexG=indexG+499
}

sum(MatrixOrderPxGxC)
sum(X)
cp_decomp <- cp(as.tensor(X), num_components = 6)
str(cp_decomp$U)


data_frame(component = c(rep("u[1]", 32),
                         rep("v[1]", 482),
                         rep("w[1]", 499)),
           value = c(cp_decomp$U[[1]][ , 1],
                     cp_decomp$U[[2]][ , 1], 
                     cp_decomp$U[[3]][ , 1] ),
           index = c((1:32), (1:482), (1:499))) %>%
  ggplot(aes(index, value)) + geom_line() +
  facet_wrap(~component, scales = "free", nrow = 3,
             labeller = labeller(component = label_parsed)) +
  theme(axis.title = element_blank())


pdf( file = paste0("temp/allComponents",k,".pdf"),  onefile = TRUE, width = 9, height = 7)
data_frame(component = c(rep("u[1]", 482), rep("u[2]", 482), rep("u[3]", 482),
                         rep("v[1]", 499), rep("v[2]", 499), rep("v[3]", 499),
                         rep("w[1]", 32), rep("w[2]", 32), rep("w[3]", 32)),
           value = c(cp_decomp$U[[1]][ , 1], cp_decomp$U[[1]][ , 2], cp_decomp$U[[1]][ , 3],
                     cp_decomp$U[[2]][ , 1], cp_decomp$U[[2]][ , 2], cp_decomp$U[[2]][ , 3],
                     cp_decomp$U[[3]][ , 1], cp_decomp$U[[3]][ , 2], cp_decomp$U[[3]][ , 3]),
           index = c(rep(1:482, 3),
                     rep(1:499, 3), 
                     rep(1:32, 3))) %>%
  ggplot(aes(index, value)) + geom_line() +
  facet_wrap(~component, scales = "free", nrow = 3,
             labeller = labeller(component = label_parsed)) +
  theme(axis.title = element_blank())
dev.off()