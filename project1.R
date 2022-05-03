source("libraries.R", encoding="UTF-8")

load("Melon.rdata")
set.seed(0875362)
mygroup <- which(rmultinom(1, 1, c(0.25,0.25,0.25,0.25)) == 1)
mysample <- sample(which(y == mygroup), 180)
X_train <- data.frame(X[mysample[1:90], ])
X_valid <- data.frame(X[mysample[91:180], ])
colorCant <- c("blue3", "orange", "darkgreen","red")

#Here we plot the spectra of the whole dataset and the training data

matplot(t(X), lty=1, type="l",col = factor(y),xlab = "Wavelength",ylab="",main = "Spectra")
legend(x=150,y=3,legend=c("D","H","Ha","E"),col=unique(factor(y)),lty=1)
matplot(t(X_train), lty=1, type="l",col = "black",xlab = "Wavelength",ylab="",main = "Training data spectra")

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

summary(train.pca)

#Both plots demonstrate the same thing - how many components to keep
plot(train.pca$eigenvalues/sum(train.pca$eigenvalues),type="b")
fviz_eig(train.pca.prcomp,choice="eigenvalue")

#select 2 components and draw some plots
train.pca.2 <- PcaClassic(X_train, k=2, scale=T)
train.pca.2.prcomp <- getPrcomp(train.pca.2)
fviz_pca_biplot(train.pca.2.prcomp, label="var")

matplot(train.pca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings", main="Classical PCA loadings", 
        ylim=c(-0.1, 0.2), yaxt="n")
axis(2, at=seq(-0.1, 0.2, 0.05), labels=seq(-0.1, 0.2, 0.05))
legend("topleft",legend = c("v1","v2"),col=c("black","red"),lty=c(1,2))

plot(train.pca.2$scores,col=y,asp=1, pch=19, main="Classic PCA scores")

#QUESTION 3

#Robust PCA
train.robpca.2 <- PcaHubert(X_train,k=2,scale = mad)
outliers.robust <- unique(c(which(train.robpca.2$od > train.robpca.2$cutoff.od),which(train.robpca.2$sd > train.robpca.2$cutoff.sd))) %>% sort()
outliers.robust.sd <- train.robpca.2$sd[outliers.robust]
outliers.robust.od <- train.robpca.2$od[outliers.robust]

matplot(train.robpca.2$loadings, type="l", xlab="Wavelength", ylab="Loadings",main="Robust PCA loadings")
legend("bottomleft",legend = c("v1","v2"),col=c("black","red"),lty=c(1,2))

plot(train.robpca.2$scores,asp=1, main="Robust PCA scores",pch=19)

#investigate two groups - very clear division with two outliers (which we detect later)
possibleGroup <-which(X_train[,100]>1)
plot(train.robpca.2$scores,pch=19,col="blue", main = "Score plot with two groups marked")
points(train.robpca.2$scores[possibleGroup,], pch=19, col="green")
legend("topleft", legend = c("group 1", "group 2"), col = c("blue", "green"), pch=19)



plot(train.robpca.2,pch=19, col = "blue", main = "Robust PCA scores with two groups detailed")
points(x = train.robpca.2$sd[possibleGroup], y = train.robpca.2$od[possibleGroup],pch=19,col="green")
legend("topright", legend = c("group 1", "group 2"), col = c("blue", "green"), pch=19)


#QUESTION 4
#validation set analysis. Here we use mainly adapted code from the last exercise session.

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

#normality tests, first of the whole cleaned training dataset and then of each of the sub-groups present
outliers.robust <- unique(c(which(train.robpca.2$od > train.robpca.2$cutoff.od),which(train.robpca.2$sd > train.robpca.2$cutoff.sd)))

X_train.scores.clean.robust <- train.robpca.2$scores[-outliers.robust,]
shapiro.test(X_train.scores.clean.robust)

plot(train.robpca.2$scores,col="blue",xlim=c(-100,60),ylim=c(-30,45),pch=19, main="Scores with classified outliers")
points(train.robpca.2$scores[outliers.robust,],col="green",pch=19)
legend("topleft",legend=c("Regular points", "Outliers"),col=c('blue',"green"),pch=19)
shapiro.test(X_train.scores.clean.robust)

#Here we perform Shapiro-Wilk tests on each of the groups found in data.

clean.in.group <- as.character(possibleGroup[is.element(as.character(possibleGroup),rownames(X_train.scores.clean.robust))])


plot(X_train.scores.clean.robust, col="blue", pch=19, main="Outlier free dataset by group")
points(X_train.scores.clean.robust[clean.in.group,],col="green", pch=19)

shapiro.test(X_train.scores.clean.robust[clean.in.group,])


shapiro.test(X_train.scores.clean.robust[!(rownames(X_train.scores.clean.robust) %in% clean.in.group),])

