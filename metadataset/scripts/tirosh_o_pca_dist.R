rm(list=ls())
author="tirosh_o"
library("FactoMineR")
library("factoextra")
load("../preprocessed/data_score_AUC_TPM_tirosh.Rdata")
load(sprintf("../preprocessed/AUCell_pca_%s.RData", author))
load(sprintf("../preprocessed/AUCell_mean_%s.RData", author))

coord <- read.table("../external_data/cell_cycle_coordinates_portal2.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)
coord <- coord[-1,]
celltype <- read.table("../external_data/cell_type_assignment_portal.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)
celltype <- celltype[-1,]
meta <- celltype[,1]
meta[meta=="0"] <- "NA"
meta[meta=="malignant"] <- "Malignant"
meta[meta=="Oligodendrocytes (non-malignant)"] <- "Oligodendrocytes"
names(meta) <- rownames(celltype)

pats <- gsub("_\\w+_\\w+$", "", rownames(celltype), perl=T)
pats <- gsub("^9", "MGH9", pats, perl=T)
names(pats) <- rownames(celltype)
allPat <- as.character(sort(unique(pats)))

colnames(score) <- gsub("^X", "", colnames(score))
ev <- apca$svd$V
##ind <- get_pca_ind(apca)
eigen <- get_eigenvalue(apca)

for (i in rownames(score)) {
  tmp <- score[i,] - min(score[i,], na.rm=T)
  score[i,] <- (tmp / max(tmp, na.rm=T)) - amean[i]
}

pcamat <- matrix(unlist(lapply(1:ncol(score), function(cell){
  unlist(lapply(1:nrow(score), function(act, cell){sum(score[,cell] * ev[,act])}, cell=cell))
})), byrow=F, nr=nrow(score))
colnames(pcamat) <- colnames(score)

thresh=2.0
## for (thresh in c(0,1,2,3,5)) {
threshidx <- which(eigen[,"variance.percent"] > thresh)

library(parallel)
distList <- list() ## malign/malign, benign/benign distances
bList <- list() ##intra-benign distances
dtList <- list()
for (p in allPat) {
  print(p)
  pcells <- names(pats)[which(pats==p)]
  pmat <- pcamat[threshidx,pcells]
  combs <- combn(colnames(pmat), 2)
  pdist <- dist(t(pmat), "euclidean")
  distcat <- unlist(mclapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    if (meta[c1] != meta[c2])
      return("Inter")
    return(meta[c1])
  }, mc.cores=8))
  names(pdist) <- distcat
  for (d in unique(distcat)) {
    if (!d %in% names(distList))
      distList[[d]] <- list()
    distList[[d]][[p]] <- pdist[which(distcat == d)]
  }

  detailcat <- unlist(mclapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    return(paste(sort(c(meta[c1], meta[c2])), collapse=":"))
  }, mc.cores=8))
  for (x in sort(unique(detailcat))) {
    if (!x %in% names(dtList)) {
      dtList[[x]] <- pdist[which(detailcat == x)]
    } else { dtList[[x]] <- c(dtList[[x]], pdist[which(detailcat == x)])}
  }
}

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
grpo <- unlist(lapply(unique(meta), function(x){which(meta==x)}))

means <- lapply(distList, function(x){unlist(lapply(x, mean, na.rm=T))})

save(distList, file=sprintf("../external_data/distList_tirosh_o_%.1f.RData", thresh))
save(dtList, file=sprintf("../external_data/dtList_tirosh_o_%.1f.RData", thresh))

## for (thresh in c(0,1,2,3,5)) {
load(sprintf("../external_data/distList_tirosh_o_%.1f.RData", thresh))
load(sprintf("../external_data/dtList_tirosh_o_%.1f.RData", thresh))

cellorder <- c("Oligodendrocytes", "Microglia/Macrophage", "Malignant", "NA")
equivgrp <- 1:length(cellorder)
names(equivgrp) <- cellorder
patgrpmat <- rbind(meta, pats)
patgrp <- 1:length(c(cellorder, allPat))
names(patgrp) <- c(cellorder, allPat)
pgm2 <- matrix(NA, nc=ncol(patgrpmat), nr=nrow(patgrpmat))
pgm2[1,] <- patgrp[patgrpmat[1,]]
pgm2[2,] <- patgrp[patgrpmat[2,]]
grpcol <- c(rainbow[c(1,3,7,8)], viridis(length(allPat)))
names(grpcol) <- names(patgrp)
patgrpo <- order(pgm2[1,], partial=pgm2[2,])

pdf(sprintf("../plots/tirosh_o_distances_%.1f.pdf", thresh))
heatmap.2(pcamat[threshidx,patgrpo], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("PC", 1:length(threshidx), sep=""), cexRow=.6, mar=c(1,3))
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,patgrpo]), col=grpcol, axes=F)
plot.new()
legend("topleft", legend=cellorder, fill=grpcol[1:length(cellorder)], horiz=F, cex=1, box.lty=0)
legend("topright", legend=allPat, fill=grpcol[(1+length(cellorder):length(grpcol))], horiz=F, cex=.75, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(Malignant=unlist(distList[["Malignant"]]),
             Oligodendrocytes=unlist(distList[["Oligodendrocytes"]]),
             `Microglia/Macrophages`=unlist(distList[["Microglia/Macrophage"]]),
             Undefined=unlist(distList[["NA"]]),
             Inter=unlist(distList[["Inter"]])),
         col=reds8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", log="")
dev.off()

##PvClust
tummat <- pcamat[,meta[colnames(pcamat)]=="Malignant"]
tumact <- score[,meta[colnames(pcamat)]=="Malignant"]
dact <- dist(tumact)
act.hc <- hclust(dact)

library(pvclust)
pvc <- pvclust(tummat, method.dist="euclidean", method.hclust="ward.D", nboot=500, parallel=as.integer(8))
save(pvc, file="tirosh_o_malign_pvclust.RData")
pvclusters=pvpick(pvc)

origsize <- unlist(lapply(pvclusters$clusters, length))
clustmin <- 0
bigclust <- pvclusters$clusters[origsize >= clustmin]

library(gplots)
clsize <- unlist(lapply(bigclust, length))
clustcells <- unlist(bigclust)
clustcellso <- intersect(colnames(tummat)[pvc$hclust$order], clustcells)
clustid <- rep(1:length(bigclust), clsize)
names(clustid) <- unlist(bigclust)
clusto <- clustid[clustcellso]
clcol <- 1:length(clsize)%%2

## for (thresh in c(0, 1, 2, 3, 5)) {
##
threshidx <- which(eigen[,"variance.percent"] > thresh)
pc.hc <- hclust(dist(tummat[threshidx,]))
layMat <- matrix(rep(1:3, c(1,length(threshidx),nrow(score))))

pdf(sprintf("plots/tirosh_o_signif_cluster_profiles_min%d_%1.f.pdf", clustmin, thresh))
layout(layMat)
par(mar=c(0,15,0.25,0.5))
image(matrix(rep(clcol,clsize[unique(clusto)]), ncol=1), col=c("black", "grey"), axes=F)
image(t(tummat[threshidx[pc.hc$order],clustcellso]), col=bluered(100), axes=F)
axis(2, at=seq(0,1, 1/(length(threshidx)-1)), label=paste("PC", 1:length(threshidx))[pc.hc$order], cex=1, las=2)
par(mar=c(1,15,0.25,0.5))
image(t(normscore[act.hc$order,clustcellso]), col=bluered(100), axes=F)
axis(2, at=seq(0,1, 1/(nrow(score)-1)), label=rownames(score)[act.hc$order], cex=1, las=2)
dev.off()
