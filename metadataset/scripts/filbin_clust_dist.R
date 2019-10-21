rm(list=ls())
author="filbin"
load("../preprocessed/data_score_AUC_TPM_filbin.Rdata")
load(sprintf("../preprocessed/AUCell_cluster_list_%s.RData", author))

meta <- read.table("../external_data/PortalK27M_Metadata.vh20180223.txt", header=T, sep="\t", stringsAsFactors=F, row.names=1)
meta <- meta[-1,]
celltype <- meta$Type
celltype[celltype=="GBM"] <- "Malignant"
pats <- meta$Sample
names(celltype) <- names(pats) <- rownames(meta)
tsne <- read.table("../external_data/PortalK27M_tSNE.vh20180223.txt", header=T, sep="\t", stringsAsFactors=F, row.names=1)
tsne <- tsne[-1,]
cellid <- rownames(meta)#tsne only for pat data
cellid <- cellid[which(celltype[cellid] != "Filter")]

colnames(meta) <- gsub("\\.", "\\-", colnames(meta))
subtypecols <- c("OPC-variable", "OC-like", "AC-like", "OPC-like")

subtype <- unlist(apply(meta[,subtypecols], 1, function(x){
  x2=x[!is.na(x)]
  if (length(x2)==0)
    return("Undetermined")
  names(x2)[which(x2==max(x2, na.rm=T))][1]
}))

allPat <- sort(unique(gsub("^(\\w+)\\-.+$", "\\1", rownames(tsne), perl=T)))
extraPat <- setdiff(as.character(sort(unique(pats[cellid]))), c("Oligo",allPat))
extraPat <- extraPat[grep("_", extraPat)]
superPat <- unique(gsub("_.+$", "", extraPat))
gbmPat <- setdiff(as.character(sort(unique(pats[cellid]))), c("Oligo",allPat,extraPat))
allincPat <- c(allPat, extraPat, gbmPat)

##ind <- get_pca_ind(metapca)
colnames(score) <- gsub("\\.", "\\-", colnames(score))

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
distList <- subList <- dtList <- dtList2 <- list()
lenmal <- malv <- subv <- c()
for (p in c(allincPat)) {
  print(p)
  pcells <- names(pats[cellid])[which(pats[cellid]==p)]
  malcells <- names(pats[cellid])[which(pats[cellid]==p & celltype[cellid] == "Malignant")]
  lenmal[p] <- length(malcells)
  pmat <- clMat[,pcells]
  combs <- combn(colnames(pmat), 2)
  pdist <- dist(t(pmat), "euclidean")
  malv <- c(malv, dist(t(pmat[,malcells]), "euclidean"))
  subv <- c(subv, dist(meta[malcells,subtypecols], method="euclid"))

  distcat <- unlist(mclapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    if (celltype[c1] != celltype[c2])
      return("Inter")
    return(celltype[c1])
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
    return(paste(sort(c(celltype[c1], celltype[c2])), collapse=":"))
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
grpo <- unlist(lapply(unique(celltype), function(x){which(celltype==x)}))

means <- lapply(distList, function(x){unlist(lapply(x, mean, na.rm=T))})

save(distList, file=sprintf("../external_data/distList_filbin_%d.RData", nbest))
save(dtList, file=sprintf("../external_data/dtList_filbin_%d.RData", nbest))

## for (nbest in 6:10) {
load(sprintf("../external_data/distList_filbin_%d.RData", nbest))
load(sprintf("../external_data/dtList_filbin_%d.RData", nbest))

cellorder <- c("Oligodendrocyte", "Immune cell", "PDX", "Cellline", "Malignant")
equivgrp <- 1:length(cellorder)
names(equivgrp) <- cellorder
patgrpmat <- rbind(celltype[cellid], pats[cellid])
patgrp <- 1:length(c(cellorder, allincPat))
names(patgrp) <- c(cellorder, allincPat)
pgm2 <- matrix(NA, nc=ncol(patgrpmat), nr=nrow(patgrpmat))
pgm2[1,] <- patgrp[patgrpmat[1,]]
pgm2[2,] <- patgrp[patgrpmat[2,]]
grpcol <- c(rainbow[c(1,3,4,5,7)], viridis(length(allincPat)))
names(grpcol) <- names(patgrp)
patgrpo <- order(pgm2[1,], partial=pgm2[2,])

pdf(sprintf("../plots/filbin_distances_%d.pdf", nbest))
heatmap.2(clMat[,patgrpo], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("CL", 1:nbest, sep=""), cexRow=.6, mar=c(1,3))
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,patgrpo]), col=grpcol, axes=F)
plot.new()
legend("topleft", legend=cellorder, fill=grpcol[1:length(cellorder)], horiz=F, cex=1, box.lty=0)
legend("topright", legend=allincPat, fill=grpcol[(1+length(cellorder):length(grpcol))], horiz=F, cex=.75, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(Malignant=unlist(distList[["Malignant"]]),
             Oligodendrocyte=unlist(distList[["Oligodendrocyte"]]),
             `Immune cell`=unlist(distList[["Immune cell"]]),
             PDX=unlist(distList[["PDX"]]),
             Cellline=unlist(distList[["Cellline"]]),
             Inter=unlist(distList[["Inter"]])),
         col=purples8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", log="")

## Special for patient BCH869
par(mar=c(6,4,.5,.5), mfrow=c(1,1))
p=superPat
pList <- list(Malignant=unlist(distList[["Malignant"]][[p]]),
              DGC100=unlist(distList[["Cellline"]][["BCH869_DGC100"]]),
              DGC75=unlist(distList[["Cellline"]][["BCH869_DGC75"]]),
              GF=unlist(distList[["Cellline"]][["BCH869_GF"]]),
              GS=unlist(distList[["Cellline"]][["BCH869_GS"]]),
              SF=unlist(distList[["Cellline"]][["BCH869_SF"]]),
              PDX=unlist(distList[["PDX"]][["BCH869_PDX"]]))
beanplot(pList, col=purples8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", log="", las=2)
dev.off()

deriv <- setdiff(names(pList), "Malignant")
derivMat <- matrix(NA, nr=2, nc=length(deriv))
colnames(derivMat) <- deriv
rownames(derivMat) <- c("t", "wilcox")
for (d in deriv) {
  derivMat[,d] <- c(t.test(pList[["Malignant"]], pList[[d]])$p.value,
                    wilcox.test(pList[["Malignant"]], pList[[d]])$p.value)
}
write.table(derivMat, file=sprintf("../external_data/BCH869_distribs_%d.txt", nbest), col.names=NA, sep="\t", quote=F)
