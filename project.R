source("libraries.R", encoding="UTF-8")

load("Melon.rdata")
set.seed(0875362)
mygroup <- which(rmultinom(1, 1, c(0.25,0.25,0.25,0.25)) == 1)
mysample <- sample(which(y == mygroup), 180)
X_train <- data.frame(X[mysample[1:90], ])
X_valid <- data.frame(X[mysample[91:180], ])
#matplot(t(X[1300:1600,]), lty=1, type="l")

#Check column variances to see if we need to scale variables:
X_train %>% as.data.frame() %>% summarise_if(is.numeric, var)

#as some variables have much larger variances, I think we need to scale them before applying PCA.

train.pca <- PcaClassic(X_train, scale = T)
train.pca.prcomp <- getPrcomp(train.pca)
summary(train.pca)
plot(train.pca$eigenvalues/sum(train.pca$eigenvalues),type="b")
fviz_eig(train.pca.prcomp,choice="eigenvalue")


#train.pca.ns <- PcaClassic(X_train)
#plot(train.pca.ns$eigenvalues/sum(train.pca.ns$eigenvalues))
#summary(train.pca.ns)

#select 2 components
train.pca.2 <- PcaClassic(X_train, k=2, scale=T)
plot(train.pca.2)
#biplot(train.pca.2, scale=0, cex=c(0.7,1))
#plot(train.pca,6,crit.pca.distances=0.99)

train.pca.2$scores
par(mfrow=c(2,1))
matplot(train.pca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings", main="classical PCA loadings")
matplot(t(X_train),type="l")
# train.pca.6 <- PcaClassic(X, k=6, scale = T)
# plot(train.pca.6,crit.pca.distances=0.99)
par(mfrow=c(1,1))
matplot(train.pca.2$scores,type='l')
pairs(train.pca.2$scores)#, col=y)
plot(train.pca.2$scores,col=y,asp=1, pch=19, main="Classic PCA scores")
plot(train.pca.2,crit.pca.distances=0.99, pch=19)

#Robust PCA
train.robpca.2 <- PcaHubert(X_train,k=2)
plot(train.robpca.2)
matplot(train.robpca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings",main="robust PCA loadings")
plot(train.robpca.2$scores,asp=1, main="Robust PCA scores",pch=19)


#validation set. Continue with classic PCA.
loadings.pca.2 <- as.matrix(train.pca.2$loadings)
##CHECK FOR SCALING THING. SLIDE 41,42

matr <- X_valid
matr.avg <- colMeans(matr)
validation.scores <- as.matrix(matr - rep(1,90) %*% t(matr.avg)) %*% loadings.pca.2
predicted_values <- validation.scores %*% t(loadings.pca.2) + rep(1,90) %*% t(matr.avg)
valid.pca.2 <- PcaClassic(X_valid,k=2,scale = T)

scores.valid <- as.matrix(valid.pca.2$scores)
plot(train.pca.2$scores,pch=19)
points(validation.scores, col="green",pch=19)

euclnorm <- function(y) sqrt(sum(y^2))
og.distance <- apply(X_valid - predicted_values,1,euclnorm)

score.distance <- sqrt(mahalanobis(validation.scores, center = FALSE, diag(valid.pca.2$eigenvalues)))

cutoff.sd <- train.pca.2$cutoff.sd
cutoff.od <- train.pca.2$cutoff.od
ods <- c(train.pca.2$od)#, og.distance)
sds <- c(train.pca.2$sd)#, score.distance)
new.distances <- cbind(score.distance,og.distance)

xmax <- max(sds, cutoff.sd,score.distance)
ymax <- max(ods, cutoff.od,og.distance)


outl.index <- (length(ods)-3) : length(ods)
colors <- rep("black", length(ods))
colors[outl.index] <- c("red", "red", "red", "orange")
plot(sds, ods, pch = 19, xlim = c(0, xmax + 1), ylim = c(0, ymax + 0.2), col = colors,
     main = "Classical PCA", xlab = "Score distance", ylab = "Orthogonal distance")
points(new.distances,col = "blue",pch=19)
abline(h = cutoff.od, col = "red", lwd = 1.5)
abline(v = cutoff.sd, col ="red", lwd = 1.5)
