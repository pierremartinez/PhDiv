rm(list=ls())
author="li"
library("FactoMineR")
library("factoextra")
load("../preprocessed/data_score_AUC_TPM_li.Rdata")
load(sprintf("../preprocessed/AUCell_pca_%s.RData", author))
load(sprintf("../preprocessed/AUCell_mean_%s.RData", author))

li <- read.table("../external_data/GSE81861_CRC_tumor_all_cells_FPKM.csv", sep=",", stringsAsFactors=F, header=T, row.names=1)
grp <- gsub("^.+__(.+)__.+$", "\\1", colnames(li))
cellid <- gsub("__.+__.+$", "", colnames(li))
names(grp) <- cellid
rm(li)
gc()
supergrp <- grp
supergrp[!grp %in% c("NA", "Epithelial")] <- "Benign"

series <- read.table("../external_data/GSE81861_series_matrix.txt", sep="\t", stringsAsFactors=F, header=F, skip=37, fill=T, quote="\"")
series[,1] <- gsub("^!", "", series[,1])
pats <- gsub("^patient id: ", "", series[which(series[,1] == "Sample_characteristics_ch1")[1],-1])
names(pats) <- series[which(series[,1] == "Sample_title"),-1]
allPat <- unique(pats)[grep("CRC", unique(pats))]

ev <- apca$svd$V
eigen <- get_eigenvalue(apca)

for (i in rownames(score)) {
  tmp <- score[i,] - min(score[i,], na.rm=T)
  score[i,] <- (tmp / max(tmp, na.rm=T)) - amean[i]
}

pcamat <- matrix(unlist(lapply(1:ncol(score), function(cell){
  unlist(lapply(1:nrow(score), function(act, cell){sum(score[,cell] * ev[,act])}, cell=cell))
})), byrow=F, nr=nrow(score))
colnames(pcamat) <- colnames(score)

thresh=2
## for (thresh in c(0,1,2,3,5)) {
threshidx <- which(eigen[,"variance.percent"] > thresh)

distList <- list()
library(parallel)
distList <- list(Malignant=list(), Benign=list(), Undefined=list(), Inter=list()) ## malign/malign, benign/benign distances
bList <- list() ##intra-benign distances
detailList <- dtList <- list()
for (p in allPat) {
  print(p)
  pcells <- intersect(names(pats)[pats==p], cellid)
  pmat <- pcamat[,pcells,drop=F]
  if (ncol(pmat) < 2)
    next
  pcombs <- combn(colnames(pmat), 2)
  pdist <- dist(t(pmat[threshidx,]), "euclidean")
  pdistcat <- unlist(mclapply(1:ncol(pcombs), function(i){
    c1=pcombs[1,i]
    c2=pcombs[2,i]
    if (supergrp[c1] != supergrp[c2])
      return("inter")
    if (supergrp[c1] == "NA") ##unresolved
      return("u")
    if (supergrp[c1] == "Benign") ##benign
      return("b")
    return("m")
  }, mc.cores=8))
  names(pdist) <- pdistcat
  distList[["Malignant"]][[p]] <- pdist[which(pdistcat == "m")]
  distList[["Benign"]][[p]] <- pdist[which(pdistcat == "b")]
  distList[["Undefined"]][[p]] <- pdist[which(pdistcat == "u")]
  distList[["Inter"]][[p]] <- pdist[which(pdistcat == "inter")]

  ## benign only
  bidx <- which(pdistcat=="b")
  pbcat <- unlist(mclapply(bidx, function(i){ ##(1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK
    c1=pcombs[1,i]
    c2=pcombs[2,i]
    if (grp[c1] != grp[c2])
      return("Inter")
    return(grp[c1])
  }, mc.cores=8))
  pbdist <- pdist[bidx]
  names(pbdist) <- pbcat
  for (x in c("Tcell","Bcell", "Macrophage", "Fibroblast","Endothelial", "MastCell", "Inter")) {##unique(bcat)) {
    bList[[x]][[p]] <- pbdist[which(pbcat == x)]
  }

  detailcat <- unlist(mclapply(1:ncol(pcombs), function(i){
    c1=pcombs[1,i]
    c2=pcombs[2,i]
    return(paste(sort(c(grp[c1], grp[c2])), collapse=":"))
  }))
  for (x in sort(unique(detailcat))) {
    if (!x %in% names(detailList))
      detailList[[x]] <- list()
    detailList[[x]][[p]] <- pdist[which(detailcat == x)]
    if (!x %in% names(dtList)) {
      dtList[[x]] <- pdist[which(detailcat == x)]
    } else { dtList[[x]] <- c(dtList[[x]], pdist[which(detailcat == x)])}
  }
}


library(beanplot)
library(gplots)
library(beeswarm)
library(viridis)
library(RColorBrewer)
greens8 <- brewer.pal(8, "Greens")
reds8 <- brewer.pal(8, "Reds")
greens8 <- brewer.pal(8, "Greens")
purples8 <- brewer.pal(8, "Purples")
oranges8 <- brewer.pal(8, "Oranges")
greys8 <- brewer.pal(8, "Greys")
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")

save(distList, file=sprintf("../external_data/distList_li_%.1f.RData", thresh))
save(bList, file=sprintf("../external_data/bList_li_%.1f.RData", thresh))
load(sprintf("../external_data/distList_li_%.1f.RData", thresh))
load(sprintf("../external_data/bList_li_%.1f.RData", thresh))

meanm <- unlist(lapply(distList[["Malignant"]], mean, na.rm=T))
meanb <- unlist(lapply(distList[["Benign"]], mean, na.rm=T))

patgrpmat <- rbind(grp[cellid],pats[cellid])
cellorder <- c("Tcell", "Bcell", "Macrophage", "Fibroblast", "Endothelial", "MastCell", "Malignant", "NA")
patgrp <- 1:length(c(cellorder, allPat))
names(patgrp) <- c(cellorder, allPat)
pgm2 <- matrix(NA, nc=ncol(patgrpmat), nr=nrow(patgrpmat))
pgm2[1,] <- patgrp[gsub("Epithelial", "Malignant", patgrpmat[1,])]
pgm2[2,] <- patgrp[patgrpmat[2,]]
grpcol <- c(rainbow, viridis(length(allPat)))
names(grpcol) <- names(patgrp)
pato <- unlist(lapply(allPat, function(x){which(pats[cellid]==x)}))
patgrpo <- unlist(lapply(cellorder, function(x){which(grp[names(pato)]==gsub("Malignant", "Epithelial", x))}))

#matrix(rep(1:length(grpcol), unlist(lapply(names(grpcol),
                                           
pdf(sprintf("../plots/li_distances_%.1f.pdf", thresh))
heatmap.2(pcamat[threshidx,patgrpo], Colv=F, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("PC", 1:length(threshidx), sep=""), mar=c(1,9), cexRow=.6)
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,patgrpo]), col=grpcol, axes=F)
plot.new()
legend("left", legend=cellorder, fill=grpcol[cellorder], horiz=F, cex=1, box.lty=0)
legend("right", legend=allPat, fill=grpcol[allPat], horiz=F, cex=1, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(Malignant=unlist(distList[["Malignant"]]),
              Tcell=unlist(bList[["Tcell"]]),
              Bcell=unlist(bList[["Bcell"]]),
              Macrophage=unlist(bList[["Macrophage"]]),
              Fibroblast=unlist(bList[["Fibroblast"]]),
              ##Endothelial=unlist(bList[["Endothelial"]]), ## Removed because n=1
              Inter=unlist(bList[["Inter"]])),
         col=greens8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", las=2, log="")
dev.off()

## Compare with Normal
nli <- read.table("../external_data/GSE81861_CRC_NM_all_cells_FPKM.csv", sep=",", stringsAsFactors=F, header=T, row.names=1)
ngrp <- gsub("^.+__(.+)__.+$", "\\1", colnames(nli))
ncellid <- gsub("__.+__.+$", "", colnames(nli))
npats <- gsub("^.+__.+__", "", colnames(nli))
names(ngrp) <- names(npats) <- ncellid
rm(nli)
gc()
nsupergrp <- ngrp
nsupergrp[!ngrp %in% c("NA", "Epithelial")] <- "Benign"

stmp = score
load("../preprocessed/data_score_AUC_TPM_li2.Rdata")
nscore=score
score=stmp

for (i in rownames(nscore)) {
  tmp <- nscore[i,] - min(nscore[i,], na.rm=T)
  nscore[i,] <- (tmp / max(tmp, na.rm=T)) - amean[i]
}

npcamat <- matrix(unlist(lapply(1:ncol(nscore), function(cell){
  unlist(lapply(1:nrow(nscore), function(act, cell){sum(nscore[,cell] * ev[,act])}, cell=cell))
})), byrow=F, nr=nrow(nscore))
colnames(npcamat) <- colnames(nscore)

ndistList <- list()
library(parallel)
ndistList <- list(Epithelial=list(), Benign=list(), Undefined=list(), Inter=list()) ## malign/malign, benign/benign distances
nbList <- list() ##intra-benign distances
ndetailList <- ndtList <- list()
for (p in allPat) {
  print(p)
  npcells <- intersect(names(pats)[pats==p], ncellid)
  pmat <- npcamat[,npcells,drop=F]
  if (ncol(pmat) < 2)
    next
  pcombs <- combn(colnames(pmat), 2)
  pdist <- dist(t(pmat[threshidx,]), "euclidean")
  pdistcat <- unlist(mclapply(1:ncol(pcombs), function(i){
    c1=pcombs[1,i]
    c2=pcombs[2,i]
    if (nsupergrp[c1] != nsupergrp[c2])
      return("inter")
    if (nsupergrp[c1] == "NA") ##unresolved
      return("u")
    if (nsupergrp[c1] == "Benign") ##benign
      return("b")
    return("e")
  }, mc.cores=8))
  names(pdist) <- pdistcat
  ndistList[["Epithelial"]][[p]] <- pdist[which(pdistcat == "e")]
  ndistList[["Benign"]][[p]] <- pdist[which(pdistcat == "b")]
  ndistList[["Undefined"]][[p]] <- pdist[which(pdistcat == "u")]
  ndistList[["Inter"]][[p]] <- pdist[which(pdistcat == "inter")]

  ## benign only
  bidx <- which(pdistcat=="b")
  pbcat <- unlist(mclapply(bidx, function(i){ ##(1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK
    c1=pcombs[1,i]
    c2=pcombs[2,i]
    if (ngrp[c1] != ngrp[c2])
      return("Inter")
    return(ngrp[c1])
  }, mc.cores=8))
  pbdist <- pdist[bidx]
  names(pbdist) <- pbcat
  for (x in c("Tcell","Bcell", "Macrophage", "Fibroblast","Endothelial", "MastCell", "Inter")) {##unique(bcat)) {
    nbList[[x]][[p]] <- pbdist[which(pbcat == x)]
  }

  ndetailcat <- unlist(mclapply(1:ncol(pcombs), function(i){
    c1=pcombs[1,i]
    c2=pcombs[2,i]
    return(paste(sort(c(ngrp[c1], ngrp[c2])), collapse=":"))
  }))
  for (x in sort(unique(ndetailcat))) {
    if (!x %in% names(ndetailList))
      ndetailList[[x]] <- list()
    ndetailList[[x]][[p]] <- pdist[which(ndetailcat == x)]
    if (!x %in% names(ndtList)) {
      ndtList[[x]] <- pdist[which(ndetailcat == x)]
    } else { ndtList[[x]] <- c(ndtList[[x]], pdist[which(ndetailcat == x)])}
  }
}

nmeane <- unlist(lapply(ndistList[["Epithelial"]], mean, na.rm=T))
nmeanb <- unlist(lapply(ndistList[["Benign"]], mean, na.rm=T))
nmeant <- unlist(lapply(nbList[["Tcell"]], mean, na.rm=T))

npatgrpmat <- rbind(ngrp[ncellid], pats[ncellid])
ncellorder <- c("Tcell", "Bcell", "Macrophage", "Fibroblast", "Endothelial", "MastCell", "Epithelial", "NA")
npatgrp <- 1:length(c(ncellorder, allPat))
names(npatgrp) <- c(ncellorder, allPat)
pgm2 <- matrix(NA, nc=ncol(npatgrpmat), nr=nrow(npatgrpmat))
pgm2[1,] <- npatgrp[npatgrpmat[1,]]
pgm2[2,] <- npatgrp[npatgrpmat[2,]]
ngrpcol <- c(rainbow, viridis(length(allPat)))
names(ngrpcol) <- names(npatgrp)
npato <- unlist(lapply(allPat, function(x){which(pats[ncellid]==x)}))
npatgrpo <- unlist(lapply(ncellorder, function(x){which(ngrp[names(npato)]==x)}))

pdf(sprintf("../plots/li_normal_distances_%.1f.pdf", thresh))
heatmap.2(npcamat[threshidx,npatgrpo], Colv=F, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("PC", 1:length(threshidx), sep=""), mar=c(1,9), cexRow=.6)
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,npatgrpo]), col=ngrpcol, axes=F)
plot.new()
legend("left", legend=ncellorder, fill=ngrpcol[ncellorder], horiz=F, cex=1, box.lty=0)
legend("right", legend=allPat, fill=ngrpcol[allPat], horiz=F, cex=1, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(Epithelial=unlist(ndistList[["Epithelial"]]),
              Tcell=unlist(nbList[["Tcell"]]),
              Bcell=unlist(nbList[["Bcell"]]),
              Macrophage=unlist(nbList[["Macrophage"]]),
              Fibroblast=unlist(nbList[["Fibroblast"]]),
              ##Endothelial=unlist(nbList[["Endothelial"]]), ## removed because n=0
              Inter=unlist(nbList[["Inter"]])),
         col=greens8[c(8,1,3,5)], cut=.1, what=c(1,1,1,0), ylab="Distance", las=2, log="")

par(mar=c(6,4,2,.5), mfrow=c(2,2))
for (p in allPat) {
  mtab <- table(grp[names(pats[cellid])[pats[cellid]==p]])
  ntab <- table(ngrp[names(pats[ncellid])[pats[ncellid]==p]])
  mtab <- mtab[mtab>1]
  ntab <- ntab[ntab>1]
  ov <- intersect(names(mtab), names(ntab))
  
  for (cat in ov) {
    if (cat == "Epithelial") {
      wt <- wilcox.test(unlist(ndistList[["Epithelial"]][[p]]), unlist(distList[["Malignant"]][[p]]))
      beanplot(list(unlist(ndistList[["Epithelial"]][[p]]), unlist(distList[["Malignant"]][[p]])),
               col=greens8[c(3,1,4,6)], cut=.1, what=c(1,1,1,0), ylab="Distance", cex.axis=1, main=p,
               sub=sprintf("p=%.3f (Wilcoxon test)", wt$p.value), names=c("Normal Epithelial", "Malignant Epithelial"), log="")
    } else {
      wt <- wilcox.test(unlist(nbList[[cat]][[p]]), unlist(bList[[cat]][[p]]))
      beanplot(list(unlist(nbList[[cat]][[p]]), unlist(bList[[cat]][[p]])),
               col=greens8[c(3,1,4,6)], cut=.1, what=c(1,1,1,0), ylab="Distance", cex.axis=1, main=p,
               sub=sprintf("p=%.3f (Wilcoxon test)", wt$p.value), names=c(paste(c("Normal", "Cancer"), cat)), log="")
    }
  }
}
dev.off()
