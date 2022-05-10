library(ggplot2)
library(tidyverse)
library(reshape2)
library(rrcov)
library(robustbase)
library(corrplot)
library(factoextra)
library(ggfortify)
library(pheatmap)
library(cluster)

load("Melon.rdata")
set.seed(0875362)
alldata <- data.frame(cbind(X,y))
mydata <- alldata %>% group_by(y) %>% sample_n(50)
X <- mydata[, -257]
group <- mydata$y
# colors for the 4 known groups
colors.gr <- c("blue","red","darkgreen","orange")
# colors for the clusters
colors.cl <- c("cadetblue","firebrick","darkolivegreen3","gold","brown","purple")

#QUESTION 1
par(mfrow=c(2,2))
matplot(t(mydata[mydata$y==1,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 1", xlab = "Wavelength", col = colors.gr[1])
matplot(t(mydata[mydata$y==2,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 2", xlab = "Wavelength", col = colors.gr[2])
matplot(t(mydata[mydata$y==3,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 3", xlab = "Wavelength", col = colors.gr[3])
matplot(t(mydata[mydata$y==4,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 4", xlab = "Wavelength", col = colors.gr[4])
#mtext("Spectral plots by group", side = 3, line = -3, outer = TRUE)

par(mfrow=c(1,1))

pheatmap(X[order(mydata$y),], legend = TRUE, cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,color=COL2("BrBG"))

#QUESTION 2
melon.pam2 <- pam(X, k=2)
melon.pam3 <- pam(X, k=3)
melon.pam4 <- pam(X, k=4)
melon.pam5 <- pam(X, k=5)
melon.pam6 <- pam(X, k=6)

whichCluster <- melon.pam2

#confusion matrix
corrplot(table(group, whichCluster$cluster), is.corr=F, method="color",
                   tl.srt=0, tl.col="black", addgrid.col="grey", addCoef.col="grey",
                   number.cex=2, cl.pos="n", col=COL1("YlGn"))
#contingency table

contingencyTable <- table(whichCluster[["clustering"]], group)
addmargins(contingencyTable)

#mosaic plot

mosaicplot(t(contingencyTable), color = colors.cl, main = "", ylab = "cluster")

#cluster plot

#NOTE: when k=4 legends for clusters and groups merge into one (it still works ok, just instead of
#two legends we get one)

f_clust_pam <- fviz_cluster(whichCluster, geom="point", ellipse.type="convex",
                            palette=colors.cl, ellipse.border.remove = T)
f_clust_pam[["layers"]][[1]][["data"]][["cluster"]] <- as.factor(group)
f_clust_pam[["layers"]][[2]][["data"]][["cluster"]] <- as.factor(whichCluster$cluster)
f_clust_pam +
  scale_colour_manual(values=colors.gr)

#plots of spectra colored by clusters

matplot(t(X), col = factor(whichCluster$cluster), type="l", lty=1, xlab = "Wavelength", 
        ylab = "", main = "Spectral plot coloured according to clustering, k = 6")
legend("topright", legend = c(unique(whichCluster[["clustering"]])), col = unique(factor(whichCluster$cluster)),
       lty = 1, cex = 0.5, horiz = T)


#heatmap according to clustering

pheatmap(X[order(whichCluster$cluster),], legend = TRUE, cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,color=COL2("BrBG"))

#silhouette values
plot(silhouette(melon.pam2), col=colors.cl[1:2])
fviz_silhouette(melon.pam2, palette=colors.cl)

plot(silhouette(melon.pam3), col=colors.cl[1:3])
fviz_silhouette(melon.pam3, palette=colors.cl)

plot(silhouette(melon.pam4), col=colors.cl[1:4])
fviz_silhouette(melon.pam4, palette=colors.cl)

plot(silhouette(melon.pam5), col=colors.cl[1:5])
fviz_silhouette(melon.pam5, palette=colors.cl)

plot(silhouette(melon.pam6), col=colors.cl[1:6])
fviz_silhouette(melon.pam6, palette=colors.cl)

#number of clusters
fviz_nbclust(X, pam)



#QUESTION 3: hierarchical clustering

melon.agnes.avg <- agnes(X)
plot(melon.agnes.avg, which.plots=2)

melon.agnes.complete <- agnes(X, method = "complete")
melon.agnes.complete.std <- agnes(X, method = "complete", stand = T)
par(mfrow=c(1,1))
plot(melon.agnes.complete, which.plots=2)
plot(melon.agnes.complete.std, which.plots=2)

melon.agnes.complete.dg <- as.dendrogram(melon.agnes.complete)
plot(melon.agnes.complete.dg, leaflab="none")
fviz_dend(melon.agnes.complete.dg, k=6, k_colors=colors.cl, show_labels = F, 
          main = "Dendrogram of complete linkage clustering with 6 clusters")

melon.agnes.avg.dg <- as.dendrogram(melon.agnes.avg)
fviz_dend(melon.agnes.avg, k=6, k_colors = colors.cl, show_labels = F,
          main = "Dendrogram of average linkage clustering with 6 clusters")


agnes.avg_2 <- cutree(melon.agnes.avg, k=2)
agnes.cmt_2 <- cutree(melon.agnes.complete, k=2)

#Hierarchical clustering at k=2 and k-medoids at k=2 give the same results.
setequal(agnes.avg_2, agnes.cmt_2)
setequal(agnes.avg_2, melon.pam2[["clustering"]])

matplot(t(X), col = factor(cutree(melon.agnes.avg, k=2)), type="l", lty=1, 
        main="Spectral plot colored according to average linkage clustering with k = 2", ylab = "", xlab = "Wavelength")

matplot(t(X), col = factor(cutree(melon.agnes.complete, k=2)), type="l", lty=1,
        main="Spectral plot colored according to complete linkage clustering with k = 2", ylab = "", xlab = "Wavelength")

matplot(t(X), col = factor(colors.cl), type="l", lty=1,
        main="Spectral plot colored according to complete linkage clustering with k = 6", ylab = "", xlab = "Wavelength")

matplot(t(X), col = factor(colors.cl), type="l", lty=1, 
        main="Spectral plot colored according to average linkage clustering with k = 6", ylab = "", xlab = "Wavelength")


table(cutree(melon.agnes.complete, k=6), group)


