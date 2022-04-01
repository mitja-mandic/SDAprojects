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

#QUESTION 2

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
outliers <- c(14,33,42,48,88)
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

#QUESTION 3

#Robust PCA
train.robpca.2 <- PcaHubert(X_train,k=2,scale = T)
plot(train.robpca.2)
matplot(train.robpca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings",main="robust PCA loadings")
plot(train.robpca.2$scores,asp=1, main="Robust PCA scores",pch=19)

#QUESTION 4
#validation set. Continue with classic PCA.

which_pca <- train.robpca.2

X_valid.standardized <- scale(X_valid,center=which_pca$center, scale=which_pca$scale)

X_valid.scores <- X_valid.standardized %*% which_pca$loadings
X_valid.fitted <- t(t(X_valid.scores %*% t(which_pca$loadings)) + which_pca$center/which_pca$scale)

X_valid.scaled <- scale(X_valid,center = F, scale = which_pca$scale)
euclnorm <- function(y) sqrt(sum(y^2))

og.distance <- apply(X_valid.scaled -X_valid.fitted, 1, euclnorm)
score.distance <- sqrt(mahalanobis(X_valid.scores, center = FALSE, diag(which_pca$eigenvalues)))

cutoff.sd <- which_pca$cutoff.sd
cutoff.od <- which_pca$cutoff.od
ods <- c(which_pca$od)#, og.distance)
sds <- c(which_pca$sd)#, score.distance)
new.distances <- cbind(score.distance,og.distance)

xmax <- max(sds, cutoff.sd,score.distance)
ymax <- max(ods, cutoff.od,og.distance)


#outl.index <- (length(ods)-3) : length(ods)
#colors <- rep("black", length(ods))
#colors[outl.index] <- c("red", "red", "red", "orange")
plot(sds, ods, pch = 19, xlim = c(0, xmax + 1), ylim = c(0, ymax + 0.2),# col = colors,
     main = "Outlier plot of predicted and training data", xlab = "Score distance", ylab = "Orthogonal distance", col = "black")
points(new.distances,col = "blue",pch=19)
legend("topright",legend = c("Predicted values","Training data values"),col=c("blue","black"),pch=19)
abline(h = cutoff.od, col = "red", lwd = 1.5)
abline(v = cutoff.sd, col ="red", lwd = 1.5)


#QUESTION 5
X_train.scores.clean <- train.pca.2$scores[-outliers,]
shapiro.test(X_train.scores.clean)


