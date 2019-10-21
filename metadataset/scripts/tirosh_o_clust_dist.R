rm(list=ls())
author="tirosh_o"
load("../preprocessed/data_score_AUC_TPM_tirosh.Rdata")
load(sprintf("../preprocessed/AUCell_cluster_list_%s.RData", author))

coord <- read.table("../external_data/cell_cycle_coordinates_portal2.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)
coord <- coord[-1,]
celltype <- read.table("../external_data/cell_type_assignment_portal.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)
celltype <- celltype[-1,]
meta <- celltype[,1] ##cell type
meta[meta=="0"] <- "NA"
meta[meta=="malignant"] <- "Malignant"
meta[meta=="Oligodendrocytes (non-malignant)"] <- "Oligodendrocytes"
names(meta) <- rownames(celltype)

pats <- gsub("_\\w+_\\w+$", "", rownames(celltype), perl=T)
pats <- gsub("^9", "MGH9", pats, perl=T)
names(pats) <- rownames(celltype)
allPat <- as.character(sort(unique(pats)))

colnames(score) <- gsub("^X", "", colnames(score))

for (i in rownames(score)) {
  tmp <- score[i,] - min(score[i,], na.rm=T)
  score[i,] <- tmp / max(tmp, na.rm=T)
}

clusterise <- function(x, cl, n) {
  unlist(lapply(1:n, function(i) {
    mean(x[names(cl)[which(cl==i)]])
  }))
}

nbest=8
## for (nbest in 6:10) {
clbest <- clList[[as.character(nbest)]]
clMat <- apply(score, 2, clusterise, cl=clbest, n=nbest)

library(parallel)
distList <- list() ## malign/malign, benign/benign distances
bList <- list() ##intra-benign distances
dtList <- list()
for (p in allPat) {
  print(p)
  pcells <- names(pats)[which(pats==p)]
  pmat <- clMat[,pcells]
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

save(distList, file=sprintf("../external_data/distList_tirosh_o_%d.RData", nbest))
save(dtList, file=sprintf("../external_data/dtList_tirosh_o_%d.RData", nbest))

## for (nbest in 6:10) {
load(sprintf("../external_data/distList_tirosh_o_%d.RData", nbest))
load(sprintf("../external_data/dtList_tirosh_o_%d.RData", nbest))

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

pdf(sprintf("../plots/tirosh_o_distances_%d.pdf", nbest))
heatmap.2(clMat[,patgrpo], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("CL", 1:nbest, sep=""), cexRow=.6, mar=c(1,3))
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
