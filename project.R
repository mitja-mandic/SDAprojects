source("libraries.R", encoding="UTF-8")

load("Melon.rdata")
set.seed(0875362)
mygroup <- which(rmultinom(1, 1, c(0.25,0.25,0.25,0.25)) == 1)
mysample <- sample(which(y == mygroup), 180)
X_train <- data.frame(X[mysample[1:90], ])
X_valid <- data.frame(X[mysample[91:180], ])
colorCant <- c("blue3", "orange", "darkgreen","red")

matplot(t(X), lty=1, type="l",col = factor(y),xlab = "Wavelength",ylab="",main = "Spectra")
legend(x=150,y=3,legend=c("D","H","Ha","E"),col=unique(factor(y)),lty=1)

#Check column variances to see if we need to scale variables:
variances <- X_train %>% as.data.frame() %>% summarise_if(is.numeric, var)
plot(as.numeric(variances),ylab="Variance",xlab="Wavelength",main="Variances of variables")

#very interestingly my sample (training and test) only includes canteloupes of the group "Ha"
table(y[mysample])

#QUESTION 2

#as some variables have much larger variances, I think we need to scale them before applying PCA.
#work with a correlation matrix
train.pca <- PcaClassic(X_train, scale = T)
train.pca.prcomp <- getPrcomp(train.pca)
plot(train.pca$scores, main = "Classic PCA Scores", pch=19, col=factor(mysample))
legend("bottomleft",legend = cultivar_levels, col = unique(factor(mysample)),pch=19)

autoplot(train.pca.prcomp, data=X_train, colour = mygroups)
summary(train.pca)

plot(train.pca$eigenvalues/sum(train.pca$eigenvalues),type="b")
fviz_eig(train.pca.prcomp,choice="eigenvalue")


#train.pca.ns <- PcaClassic(X_train)
#plot(train.pca.ns$eigenvalues/sum(train.pca.ns$eigenvalues))
#summary(train.pca.ns)

#select 2 components
train.pca.2 <- PcaClassic(X_train, k=2, scale=T)
train.pca.2.prcomp <- getPrcomp(train.pca.2)
fviz_pca_biplot(train.pca.2.prcomp, label="var")

outliers <- c(14,33,42,48,88)
#biplot(train.pca.2, scale=0, cex=c(0.7,1))
#plot(train.pca,6,crit.pca.distances=0.99)


train.pca.2$scores

matplot(train.pca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings", main="Classical PCA loadings", 
        ylim=c(-0.1, 0.2), yaxt="n")
axis(2, at=seq(-0.1, 0.2, 0.05), labels=seq(-0.1, 0.2, 0.05))
legend("topleft",legend = c("v1","v2"),col=c("black","red"),lty=c(1,2))

matplot(t(X_train),type="l")
# train.pca.6 <- PcaClassic(X, k=6, scale = T)
# plot(train.pca.6,crit.pca.distances=0.99)
par(mfrow=c(1,1))
# matplot(train.pca.2$scores,type='l')
# pairs(train.pca.2$scores)#, col=y)
plot(train.pca.2$scores,col=y,asp=1, pch=19, main="Classic PCA scores")
#plot(train.pca.2,crit.pca.distances=0.99, pch=19)

#QUESTION 3

#Robust PCA
train.robpca.2 <- PcaHubert(X_train,k=2,scale = mad)
outliers.robust <- unique(c(which(train.robpca.2$od > train.robpca.2$cutoff.od),which(train.robpca.2$sd > train.robpca.2$cutoff.sd))) %>% sort()
outliers.robust.sd <- train.robpca.2$sd[outliers.robust]
outliers.robust.od <- train.robpca.2$od[outliers.robust]


plot(train.robpca.2,pch=19)
#text(x=train.robpca.2$sd,y=train.robpca.2$od,labels = "")
#text(x=outliers.robust.sd+0.2, y=outliers.robust.od,labels=as.character(outliers.robust))
matplot(train.robpca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings",main="Robust PCA loadings")
legend("bottomleft",legend = c("v1","v2"),col=c("black","red"),lty=c(1,2))

plot(train.robpca.2$scores,asp=1, main="Robust PCA scores",pch=19)

#QUESTION 4
#validation set.

which_pca <- train.robpca.2

X_valid.standardized <- scale(X_valid,center=which_pca$center, scale=which_pca$scale)
Dmat_X <- diag(which_pca$scale)


X_valid.scores <- as.matrix(X_valid.standardized) %*% which_pca$loadings
X_valid.fitted <- t(t(X_valid.scores %*% t(which_pca$loadings)) + which_pca$center/which_pca$scale)

plot(X_valid.scores, pch=19, main = "Predicted scores")
X_valid.scaled <- scale(X_valid,center = F, scale = which_pca$scale)
euclnorm <- function(y) sqrt(sum(y^2))

og.distance <- apply(X_valid.scaled - X_valid.fitted, 1, euclnorm)
score.distance <- sqrt(mahalanobis(X_valid.scores, center = FALSE, diag(which_pca$eigenvalues)))

cutoff.sd <- which_pca$cutoff.sd
cutoff.od <- which_pca$cutoff.od
ods <- c(which_pca$od)#, og.distance)
sds <- c(which_pca$sd)#, score.distance)
new.distances <- cbind(score.distance,og.distance)

xmax <- max(sds, cutoff.sd,score.distance)
ymax <- max(ods, cutoff.od,og.distance)

valid.robpca.2 <- PcaHubert(X_valid,k=2,scale=mad)

plot(sds, ods, pch = 19, xlim = c(0, xmax + 1), ylim = c(0, ymax + 0.2),# col = colors,
     main = "Outlier plot of predicted and training data", xlab = "Score distance", ylab = "Orthogonal distance", col = "black")
points(new.distances,col = "blue",pch=19)
#points(valid.robpca.2$sd,valid.robpca.2$od,col="green",pch=19)
legend("topright",legend = c("Predicted values","Training data values"),col=c("blue","black"),pch=19)
abline(h = cutoff.od, col = "red", lwd = 1.5)
abline(v = cutoff.sd, col ="red", lwd = 1.5)

Xhat <- X_valid.fitted %*% Dmat_X
matplot(t(Xhat),lty=1,type="l",xlab = "Wavelength",ylab="",main = "Predicted spectra, k=2",ylim = c(-2,2))
matplot(t(X_valid),type="l",lty=1, xlab = "Wavelength",ylab="", main="Observed validation data spectra",ylim = c(-2,2))



#QUESTION 5
outliers.robust <- unique(c(which(train.robpca.2$od > train.robpca.2$cutoff.od),which(train.robpca.2$sd > train.robpca.2$cutoff.sd) ))
outliers.classic <- unique(c(which(train.pca.2$od > train.pca.2$cutoff.od),which(train.pca.2$sd > train.pca.2$cutoff.sd) ))

X_train.scores.clean.classic <- train.pca.2$scores[-outliers.classic,]
shapiro.test(X_train.scores.clean.classic)

X_train.scores.clean.robust <- train.robpca.2$scores[-outliers.robust,]


plot(train.robpca.2$scores,col="blue",xlim=c(-100,60),ylim=c(-30,45),pch=19, main="Scores with classified outliers")
points(train.robpca.2$scores[outliers.robust,],col="green",pch=19)
legend("topleft",legend=c("Regular points", "Outliers"),col=c('blue',"green"),pch=19)
shapiro.test(X_train.scores.clean.robust)

