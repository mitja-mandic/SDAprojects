source("libraries.R", encoding = "UTF-8")
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

par(mfrow=c(2,2))
matplot(t(mydata[mydata$y==1,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 1", xlab = "Wavelength", col = colors.gr[1])
matplot(t(mydata[mydata$y==2,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 2", xlab = "Wavelength", col = colors.gr[2])
matplot(t(mydata[mydata$y==3,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 3", xlab = "Wavelength", col = colors.gr[3])
matplot(t(mydata[mydata$y==4,]), lty=1, type="l", ylim = c(-2,4), ylab = "", main = "group 4", xlab = "Wavelength", col = colors.gr[4])
#mtext("Spectral plots by group", side = 3, line = -3, outer = TRUE)

# x<-1:10
# par(mar=c(2.5,2.5,1,1))
# layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,3,1,3))
# plot.new()
# text(0.5,0.5,"First title",cex=2,font=2)
# plot(x)
# plot.new()
# text(0.5,0.5,"Second title",cex=2,font=2)
# hist(x)
# boxplot(x)
# barplot(x)





par(mfrow=c(1,1))

pheatmap(X[order(mydata$y),], legend = TRUE, cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,color=COL2("BrBG"))

melon.pam2 <- pam(X, k=2)
melon.pam3 <- pam(X, k=3)
melon.pam4 <- pam(X, k=4)
melon.pam5 <- pam(X, k=5)
melon.pam6 <- pam(X, k=6)

whichCluster <- melon.pam3

#confusion matrix
corrplot(table(group, whichCluster$cluster), is.corr=F, method="color",
                   tl.srt=0, tl.col="black", addgrid.col="grey", addCoef.col="grey",
                   number.cex=2, cl.pos="n", col=COL1("YlGn"))


#cluster plot
f_clust_pam <- fviz_cluster(whichCluster, geom="point", ellipse.type="convex",
                            palette=colors.cl, ellipse.border.remove = T)
f_clust_pam[["layers"]][[1]][["data"]][["cluster"]] <- as.factor(group)
f_clust_pam[["layers"]][[2]][["data"]][["cluster"]] <- as.factor(whichCluster$cluster)
f_clust_pam +
  scale_colour_manual(values=colors.gr)

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

#plots of spectra colored by clusters

matplot(t(X), col = factor(whichCluster$cluster), type="l", lty=1)

#heatmap according to clustering

pheatmap(X[order(whichCluster$cluster),], legend = TRUE, cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,color=COL2("BrBG"))


melon.agnes.avg <- agnes(X)
plot(melon.agnes.avg, which.plots=2)

melon.agnes.complete <- agnes(X, method = "complete")
melon.agnes.complete.std <- agnes(X, method = "complete", stand = T)
par(mfrow=c(1,1))
plot(melon.agnes.complete, which.plots=2)
plot(melon.agnes.complete.std, which.plots=2)

melon.agnes.complete.dg <- as.dendrogram(melon.agnes.complete)
plot(melon.agnes.complete.dg, leaflab="none")

melon.agnes.complete.std.dg <- as.dendrogram(melon.agnes.complete.std)
plot(melon.agnes.complete.std.dg, leaflab="none")

cutree(melon.agnes.avg, k=2)
cutree(melon.agnes.complete, k=2)

matplot(t(X), col = factor(cutree(melon.agnes.avg, k=2)), type="l", lty=1, main="Average")
matplot(t(X), col = factor(cutree(melon.agnes.complete, k=2)), type="l", lty=1, main = "Complete", palette(colors.cl))



