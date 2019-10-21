rm(list=ls())
author="venteicher"
load("../preprocessed/data_score_AUC_TPM_venteicher.Rdata")
load(sprintf("../preprocessed/AUCell_cluster_list_%s.RData", author))

meta <- read.table("../external_data/IDH_A_cell_type_assignment_portal_v2.txt", header=T, sep="\t", stringsAsFactors=F, row.names=1)
meta <- meta[-1,]
celltype <- meta$CLUSTER
celltype[celltype=="malignant"] <- "Malignant"
celltype[celltype=="microglia/macrophage"] <- "Microglia/Macrophage"
celltype[celltype=="oligodendrocytes"] <- "Oligodendrocytes"
celltype[celltype=="T cells"] <- "Tcells"
pats <- meta$Tumor.index
names(celltype) <- names(pats) <-  gsub(" ", "", rownames(meta))
cellid <- gsub(" ", "", rownames(meta))

allPat <- sort(unique(pats))

colnames(score) <- gsub("\\.", "\\-", colnames(score))
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
clMat <- clMat[,cellid]

library(parallel)
distList <- list() ## malign/malign, benign/benign distances
bList <- list() ##intra-benign distances
dtList <- list()
for (p in c(allPat)) {
  print(p)
  pcells <- names(pats[cellid])[which(pats[cellid]==p)]
  pmat <- clMat[,pcells]
  combs <- combn(colnames(pmat), 2)
  pdist <- dist(t(pmat), "euclidean")
  distcat <- unlist(lapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    if (celltype[c1] != celltype[c2])
      return("Inter")
    return(celltype[c1])
  }))
  names(pdist) <- distcat
  for (d in unique(distcat)) {
    if (!d %in% names(distList))
      distList[[d]] <- list()
    distList[[d]][[p]] <- pdist[which(distcat == d)]
  }

  detailcat <- unlist(lapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    return(paste(sort(c(celltype[c1], celltype[c2])), collapse=":"))
  }))
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
grpo <- unlist(lapply(unique(celltype), function(x){which(celltype==x)}))

means <- lapply(distList, function(x){unlist(lapply(x, mean, na.rm=T))})

save(distList, file=sprintf("../external_data/distList_venteicher_%d.RData", nbest))
save(dtList, file=sprintf("../external_data/dtList_venteicher_%d.RData", nbest))

## for (nbest in 6:10) {
load(sprintf("../external_data/distList_venteicher_%d.RData", nbest))
load(sprintf("../external_data/dtList_venteicher_%d.RData", nbest))

cellorder <- c("Oligodendrocytes", "Tcells", "Microglia/Macrophage", "Malignant")
equivgrp <- 1:length(cellorder)
names(equivgrp) <- cellorder
patgrpmat <- rbind(celltype[cellid], pats[cellid])
patgrp <- 1:length(c(cellorder, allPat))
names(patgrp) <- c(cellorder, allPat)
pgm2 <- matrix(NA, nc=ncol(patgrpmat), nr=nrow(patgrpmat))
pgm2[1,] <- patgrp[patgrpmat[1,]]
pgm2[2,] <- patgrp[patgrpmat[2,]]
grpcol <- c(rainbow[c(1,2,3,7)], viridis(length(allPat)))
names(grpcol) <- names(patgrp)
patgrpo <- order(pgm2[1,], partial=pgm2[2,])

meanm <- lapply(distList[["Malignant"]], mean)

pdf(sprintf("../plots/venteicher_distances_%d.pdf", nbest))
heatmap.2(clMat[,patgrpo], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("CL", 1:nbest, sep=""), cexRow=.6, mar=c(1,3))
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,patgrpo]), col=grpcol, axes=F)
plot.new()
legend("topleft", legend=cellorder, fill=grpcol[1:length(cellorder)], horiz=F, cex=1, box.lty=0)
legend("topright", legend=allPat, fill=grpcol[(1+length(cellorder):length(grpcol))], horiz=F, cex=.75, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(Malignant=unlist(distList[["Malignant"]]),
              Tcells=unlist(distList[["Tcells"]]),
              Oligodendrocytes=unlist(distList[["Oligodendrocytes"]]),
              `Microglia/Macrophage`=unlist(distList[["Microglia/Macrophage"]]),
              Inter=unlist(distList[["Inter"]])),
         col=greys8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", log="")
dev.off()
