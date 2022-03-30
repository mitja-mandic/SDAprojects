source("libraries.R", encoding="UTF-8")

load("Melon.rdata")
set.seed(0875362)
mygroup <- which(rmultinom(1, 1, c(0.25,0.25,0.25,0.25)) == 1)
mysample <- sample(which(y == mygroup), 180)
X_train <- data.frame(X[mysample[1:90], ])
X_valid <- data.frame(X[mysample[91:180], ])
matplot(t(X[1300:1600,]), lty=1, type="l")

train.pca <- PcaClassic(X, scale = T)
summary(train.pca)
plot(train.pca$eigenvalues/sum(train.pca$eigenvalues))

train.pca.ns <- PcaClassic(X)
plot(train.pca.ns$eigenvalues/sum(train.pca.ns$eigenvalues))
summary(train.pca.ns)

#select 6 components
train.pca.ns.6 <- PcaClassic(X, k=6)
biplot(train.pca.ns.6, scale=0)#, cex=c(0.7,1))

train.pca.ns.6$scores
matplot(train.pca.ns.6$loadings, type="l", xlab="Wavelength", ylab="Loadings")
