setwd("../preprocessed")
mats <- dir("./", pattern="data.score.*Rdata")
sets <- gsub("data_score_AUC_TPM_(\\w+).Rdata", "\\1", mats, perl=T)
setnames <- c("fan", "filbin", "li_tumour", "li_normal", "neftel", "patel", "tirosh_m", "tirosh_o", "venteicher")
authors <- c("fan", "filbin", "li", "li", "neftel", "patel", "tirosh_m", "tirosh_o", "venteicher")

library(beanplot)
library(beeswarm)
library(viridis)
library(RColorBrewer)
library(gplots)
blues8 <- brewer.pal(8, "Blues")
reds8 <- brewer.pal(8, "Reds")
greens8 <- brewer.pal(8, "Greens")
purples8 <- brewer.pal(8, "Purples")
oranges8 <- brewer.pal(8, "Oranges")
greys8 <- brewer.pal(8, "Greys")

rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")

setcol <- c("gold", purples8[5], greens8[5], greens8[3], oranges8[5], "black", blues8[5], reds8[5], greys8[5])

metaMat <- lens <- c()
for (m in mats) {
  mat <- get(load(m))
  lens <- c(lens, ncol(mat))
  for (i in 1:nrow(mat)) {
    tmp <- mat[i,] - min(mat[i,], na.rm=T)
    mat[i,] <- tmp / max(tmp, na.rm=T)
  }
  metaMat <- cbind(metaMat, mat)
}
save(metaMat, file="AUCell_metaMat.RData")
setv <- rep(sets, lens)
authv <- rep(authors, lens)

corMat <- matrix(NA, nr=nrow(metaMat), nc=nrow(metaMat))
colnames(corMat) <- rownames(corMat) <- rownames(metaMat)
for (i in 1:(nrow(metaMat)-1)) {
  for (j in (i+1):nrow(metaMat)) {
    cor <- cor.test(metaMat[i,], metaMat[j,])$estimate
    corMat[i,j] <- corMat[j,i] <- cor
   }
}
library(gplots)
pdf("../plots/metaMat_corrMat.pdf")
heatmap.2(corMat, trace="n", col=bluered(100), mar=c(9,9), cexRow=.6, cexCol=.6)
dev.off()

pdf("../plots/metaMat_metaMat.pdf")
metahm <- heatmap.2(metaMat, trace="n", col=bluered(100), mar=c(1,9), cexRow=.6, cexCol=.6, labCol=NA, plot=F)
par(mfrow=c(3,1))
image(matrix(as.numeric(as.factor(setv[order(order(metahm$colInd))])), ncol=1), col=setcol, axes=F)
image(matrix(as.numeric(as.factor(setv[])), ncol=1), col=setcol, axes=F)
plot.new()
legend("center", legend=setnames, fill=setcol, horiz=T, cex=0.7, bty="n")
dev.off()

layMat <- matrix(c(rep(rep(c(1,2,3), c(1,5,2)), 1),
                   rep(rep(c(4,5,6), c(1,5,2)), 8),
                   rep(rep(c(7,8,9), c(1,5,2)), 1),
                   rep(10,8)), byrow=T, nc=8)
png("AUCell_metaMat.png", res=1200, pointsize=2.5)
layout(layMat)
par(mar=c(.5,.25,.5,.25))
image(t(matrix(1:100, nr=1)), col=bluered(100), yaxt="n", xaxt="n", xlab="Normalised activity score", bty="n")
axis(1, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=0.3, tick=F, pos=3.75)
par(mar=c(0,0,0,0))
plot(metahm$colDendrogram, leaflab="none", yaxt="n", edgePar=list(lwd=0.025))
plot.new()
plot(metahm$rowDendrogram, leaflab="none", yaxt="n", horiz=T, edgePar=list(lwd=0.1))
par(mar=c(0.45, 0.35, 0.45, 0.35))
image(t(metaMat)[order(order(metahm$colInd)),], col=bluered(100), axes=F)
par(mar=c(0.45, 0.05, 0.45, 0.05))
plot(c(NA),c(NA), ylim=c(1,nrow(metaMat)), xlim=c(0,1), axes=F)
text(x=rep(0,nrow(metaMat)), y=(1:nrow(metaMat)-1)*1.055, labels=rownames(metaMat), adj=0, cex=0.3)
par(mar=c(0.45, 0.35, 0.45, 0.35))
plot.new()
image(matrix(as.numeric(as.factor(setv[order(order(metahm$colInd))])), ncol=1), col=setcol, axes=F)
plot.new()
par(mar=c(0.1,0.1, 0.1, 0.1))
plot.new()
legend("center", legend=setnames, fill=setcol, horiz=T, cex=0.3, bty="n", border=NA)
dev.off()

####################################################
########## Leave-one out clusters and PCA ##########
####################################################

library(cluster)
library(gplots)
library("ConsensusClusterPlus")
library("FactoMineR")
library("factoextra")

for (a in unique(authors)) {
  ## Clusters
  aMat <- metaMat[,which(authv != a)]
  ccp = ConsensusClusterPlus(t(aMat), maxK=15, reps=100, pItem=0.8, pFeature=.9, verbose=T,
                             title=a, clusterAlg="kmdist", distance="spearman", seed=1984.0607, plot="pdf")
  clList <- list()
  for (nbest in 2:15) {
    cl <- lapply(ccp[[nbest]]$clrs[[3]], function(x){rownames(aMat)[which(ccp[[nbest]]$clrs[[1]] == x)]})
    clv <- unlist(lapply(rownames(aMat), function(x){which(unlist(lapply(cl, function(y, x){x %in% y}, x=x)) == T)}))
    names(clv) <- rownames(aMat)
    clList[[as.character(nbest)]] <- clv
  }
  save(clList, file=sprintf("AUCell_cluster_list_%s.RData", a))

  ## PCA
  apca <- PCA(t(aMat), scale.unit = F, ncp=50, graph=F)
  var <- get_pca_var(apca)
  ind <- get_pca_ind(apca)
  eigen <- get_eigenvalue(apca)
  ev <- apca$svd$V

  amean <- apply(aMat, 1, mean, na.rm=T)
  normA <- aMat
  for (i in 1:nrow(aMat)) {
    normA[i,] <- normA[i,] - amean[i]
  }
  ##sum(ev[,1] * normA[,1])

  save(apca, file=sprintf("AUCell_pca_%s.RData", a))
  save(amean, file=sprintf("AUCell_mean_%s.RData", a))
  write.table(eigen, file=sprintf("AUCell_eigen_%s.txt", a), sep="\t", row.names=F, quote=F)

  library("corrplot")
  pdf(sprintf("AUCell_pca_%s.pdf", a))
  fviz_pca_var(apca, col.var = "black")
  fviz_pca_var(apca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE #evite le chevauchement de texte
               )
  corrplot(var$cos2, is.corr=FALSE, title="Cos2", tl.cex=0.66)
  corrplot(var$contrib, is.corr=FALSE, title="Contribution", tl.cex=0.66)
  dev.off()
}

#######################################
####    Clusters whole metaMat     ####
#######################################
library(cluster)
library(gplots)
library("ConsensusClusterPlus")
ccpmeta = ConsensusClusterPlus(t(metaMat), maxK=15, reps=100, pItem=0.8, pFeature=.9, verbose=T, 
                           title="meta", clusterAlg="kmdist", distance="spearman", seed=1984.0607, plot="pdf")

clusterise <- function(x, cl, n) {
  unlist(lapply(1:n, function(i) {
    mean(x[names(cl)[which(cl==i)]])
  }))
}

library(gplots)
pdf("../plots/metaMat_cluster_heatmap.pdf")
for (nbest in 6:10) {
  print(nbest)
  clbest <- clList[[as.character(nbest)]]
  score <- metaMat
  
  for (i in 1:nrow(score)) {
    tmp <- score[i,] - min(score[i,], na.rm=T)
    score[i,] <- tmp / max(tmp, na.rm=T)
  }
  
  clMat <- apply(score, 2, clusterise, cl=clbest, n=nbest)
  rownames(clMat) <- paste("Cluster", 1:nbest)
  
  ncorMat <- matrix(NA, nr=nrow(clMat), nc=nrow(clMat))
  colnames(ncorMat) <- rownames(ncorMat) <- rownames(clMat)
  for (i in 1:(nrow(clMat)-1)) {
    for (j in (i+1):nrow(clMat)) {
      cor <- cor.test(clMat[i,], clMat[j,], method="spearman")$estimate
      ncorMat[i,j] <- ncorMat[j,i] <- cor
    }
  }
  
  heatmap.2(clMat, trace="n", col=heat.colors(100), mar=c(1,9), cexRow=.6, cexCol=.6, labCol=NA, key.title="Cluster score")
  ##colsep=lens[1:(length(lens)-1)], sepcolor="black")
  
  heatmap.2(ncorMat, trace="n", col=bluered(100), mar=c(1,9), cexRow=1, cexCol=1, labCol=NA, key.title="Spearman's rho")
  ##colsep=lens[1:(length(lens)-1)], sepcolor="black")
  
}
dev.off()

#######################################
####       PCA whole metaMat       ####
#######################################
library("FactoMineR")
library("factoextra")
library("corrplot")
pca <- PCA(t(metaMat), scale.unit = F, ncp=50, graph=F)
var <- get_pca_var(pca)
ind <- get_pca_ind(pca)
eigen <- get_eigenvalue(pca)
ev <- pca$svd$V

mean <- apply(metaMat, 1, mean, na.rm=T)
norm <- metaMat
for (i in 1:nrow(metaMat)) {
  norm[i,] <- norm[i,] - mean[i]
}

save(pca, file="AUCell_pca_all.RData")
save(mean, file="AUCell_mean_all.RData")
write.table(eigen, file="AUCell_eigen_all.txt", sep="\t", row.names=F, quote=F)

pdf("../plots/metaMat_pca_all.pdf")
fviz_pca_var(pca, col.var = "black")
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
corrplot(var$cos2, is.corr=FALSE, title="Cos2", tl.cex=0.66)
corrplot(var$contrib, is.corr=FALSE, title="Contribution", tl.cex=0.66)
for (thresh in c(0,1,2,3,5)) {
  threshidx <- which(eigen[,"variance.percent"] > thresh)
  pcamat <- matrix(unlist(lapply(1:ncol(score), function(cell){
    unlist(lapply(1:nrow(score), function(act, cell){sum(score[,cell] * ev[,act])}, cell=cell))
  })), byrow=F, nr=nrow(score))
  colnames(pcamat) <- colnames(score)
  pcamat <- pcamat[threshidx,]
  
  pcorMat <- matrix(NA, nr=nrow(pcamat), nc=nrow(pcamat))
  colnames(pcorMat) <- rownames(pcorMat) <- rownames(pcamat)
  for (i in 1:(nrow(pcamat)-1)) {
    for (j in (i+1):nrow(pcamat)) {
      cor <- cor.test(pcamat[i,], pcamat[j,], method="spearman")$estimate
      pcorMat[i,j] <- pcorMat[j,i] <- cor
    }
  }
  
  heatmap.2(pcamat, trace="n", col=heat.colors(100), mar=c(1,9), cexRow=.6, cexCol=.6, labCol=NA, key.title="Component score")
  ##colsep=lens[1:(length(lens)-1)], sepcolor="black")
  
  heatmap.2(pcorMat, trace="n", col=bluered(100), mar=c(1,9), cexRow=1, cexCol=1, labCol=NA, key.title="Spearman's rho")
  ##colsep=lens[1:(length(lens)-1)], sepcolor="black")
}
dev.off()
