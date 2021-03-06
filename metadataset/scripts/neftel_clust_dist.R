rm(list=ls())
author="neftel"
load("../preprocessed/data_score_AUC_TPM_neftel.Rdata")
load(sprintf("../preprocessed/AUCell_cluster_list_%s.RData", author))

meta <- read.table("../external_data/IDHwtGBM.Metadata.SS2.txt", header=T, sep="\t", stringsAsFactors=F, check.names=F, row.names=1)
meta <- meta[-1,]
celltype <- meta$Cell_Type
pats <- meta$Sample
names(celltype) <- names(pats) <- rownames(meta)
cellid <- rownames(meta)
subtypecols <- c("AC-like", "MES1-like", "MES2-like", "NPC1-like", "NPC2-like", "OPC-like")

subtype <- unlist(apply(meta[,subtypecols], 1, function(x){
  x2=x[!is.na(x)]
  if (length(x2)==0)
    return("Undetermined")
  names(x2)[which(x2==max(x2, na.rm=T))][1]
}))

allPat <- sort(unique(pats))
library(gdata)
  
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
distList <- subList <- dtList <- dtList2 <- list()
lenmal <- malv <- subv <- c()
for (p in c(allPat)) {
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
  distname <- unlist(mclapply(1:ncol(combs), function(i){
      paste(c(p,combs[,i]), collapse=":")
  }, mc.cores=8))
  names(pdist) <- distname

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

  subcat <- unlist(mclapply(which(detailcat == "Malignant:Malignant"), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    if (subtype[c1] != subtype[c2])
      return("Inter")
    return(subtype[c1])
  }, mc.cores=8))

  for (d in unique(subcat)) {
    if (!d %in% names(subList))
      subList[[d]] <- list()
    subList[[d]][[p]] <- pdist[which(subcat == d)]
  }

  dt2cat <- unlist(mclapply(which(detailcat == "Malignant:Malignant"), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    return(paste(sort(c(subtype[c1], subtype[c2])), collapse=":"))
  }, mc.cores=8))
  for (x in sort(unique(dt2cat))) {
    if (!x %in% names(dtList2)) {
      dtList2[[x]] <- pdist[which(dt2cat == x)]
    } else { dtList2[[x]] <- c(dtList2[[x]], pdist[which(dt2cat == x)])}
  }
}

library(beanplot)
library(beeswarm)
library(viridis)
library(RColorBrewer)
library(gplots)
oranges8 <- brewer.pal(8, "Oranges")
reds8 <- brewer.pal(8, "Reds")
greens8 <- brewer.pal(8, "Greens")
purples8 <- brewer.pal(8, "Purples")
oranges8 <- brewer.pal(8, "Oranges")
oranges8 <- brewer.pal(8, "Oranges")
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
rainbow11 <- c("#FBB735", "#E98931", "#EB403B", "#B32E37", "#6C2A6A", "#5C4399", "#274389", "#1F5EA8", "#227FB0", "#2AB0C5", "#39C0B3")
grpo <- unlist(lapply(unique(celltype), function(x){which(celltype==x)}))

means <- lapply(distList, function(x){unlist(lapply(x, mean, na.rm=T))})

save(distList, file=sprintf("../external_data/distList_neftel_%d.RData", nbest))
save(dtList, file=sprintf("../external_data/dtList_neftel_%d.RData", nbest))
save(dtList2, file=sprintf("../external_data/dtList2_neftel_%d.RData", nbest))

## for (nbest in 6:10) {
load(sprintf("../external_data/distList_neftel_%d.RData", nbest))
load(sprintf("../external_data/dtList_neftel_%d.RData", nbest))

cellorder <- c("Oligodendrocyte", "T-cell", "Macrophage", "Malignant")
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

meanm <- unlist(lapply(distList[["Malignant"]], mean))

pdf(sprintf("../plots/neftel_distances_%d.pdf", nbest))
heatmap.2(clMat[,patgrpo], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("CL", 1:nbest, sep=""), cexRow=.6, mar=c(1,3))
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,patgrpo]), col=grpcol, axes=F)
plot.new()
legend("topleft", legend=cellorder, fill=grpcol[1:length(cellorder)], horiz=F, cex=1, box.lty=0)
legend("topright", legend=allPat, fill=grpcol[(1+length(cellorder):length(grpcol))], horiz=F, cex=.75, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(Malignant=unlist(distList[["Malignant"]]),
             Oligodendrocyte=unlist(distList[["Oligodendrocyte"]]),
             `T cell`=unlist(distList[["T-cell"]]),
             Macrophage=unlist(distList[["Macrophage"]]),
             Inter=unlist(distList[["Inter"]])),
         col=oranges8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", log="")
dev.off()

##############
## Subtypes ##
##############

smeans <- lapply(subList, function(x){unlist(lapply(x, mean, na.rm=T))})

cellorder2 <- c("AC-like", "MES1-like", "MES2-like", "NPC1-like", "NPC2-like", "OPC-like", "Undetermined")
equivgrp2 <- 1:length(cellorder2)
names(equivgrp2) <- cellorder2
patgrp2mat <- rbind(subtype[cellid], pats[cellid])
patgrp2 <- 1:length(c(cellorder2, allPat))
names(patgrp2) <- c(cellorder2, allPat)
pgm22 <- matrix(NA, nc=ncol(patgrp2mat), nr=nrow(patgrp2mat))
pgm22[1,] <- patgrp2[patgrp2mat[1,]]
pgm22[2,] <- patgrp2[patgrp2mat[2,]]
grpcol2 <- c(rainbow[c(1:length(cellorder2))], viridis(length(allPat)))
names(grpcol2) <- names(patgrp2)
patgrp2o <- order(pgm22[1,], partial=pgm22[2,])

library(vegan)
simp <- unlist(lapply(allPat, function(p){
  diversity(table(subtype[which(pats == p & celltype == "Malignant")]), index="simpson")
}))
shan <- unlist(lapply(allPat, function(p){
  diversity(table(subtype[which(pats == p & celltype == "Malignant")]))
}))

subpat <- names(lenmal)[lenmal>0]
subpatcol <- magma(length(subpat))

pdf(sprintf("../plots/neftel_subtype_distances_%d.pdf", nbest))
heatmap.2(clMat[,patgrp2o], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("CL", 1:nbest, sep=""), cexRow=.6, mar=c(1,3))
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm22[2:1,patgrp2o]), col=grpcol2, axes=F)
plot.new()
legend("topleft", legend=cellorder2, fill=grpcol2[1:length(cellorder2)], horiz=F, cex=1, box.lty=0)
legend("topright", legend=allPat, fill=grpcol2[(1+length(cellorder2):length(grpcol2))], horiz=F, cex=.75, box.lty=0)

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
ct <- cor.test(simp, meanm, method="spearman")
plot(simp, meanm, ylab="Mean malignant-malignant distance", xlab="Simpson diversity (Malignant compartment)", col=oranges8[5], pch=16,
     sub=sprintf("tau=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
abline(lm(meanm ~ simp))
ct <- cor.test(shan, meanm, method="spearman")
plot(shan, meanm, ylab="Mean malignant-malignant distance", xlab="Shannon diversity (Malignant compartment)", col=oranges8[5], pch=16,
     sub=sprintf("tau=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
abline(lm(meanm ~ shan))

par(mar=c(6,4,.5,.5), mfrow=c(1,1))
beanplot(list(`AC-like`=unlist(subList[["AC-like"]]),
             `MES1-like`=unlist(subList[["MES1-like"]]),
             `MES2-like`=unlist(subList[["MES2-like"]]),
             `NPC1-like`=unlist(subList[["NPC1-like"]]),
             `NPC2-like`=unlist(subList[["NPC2-like"]]),
             `OPC-like`=unlist(subList[["OPC-like"]]),
             ##`Undetermined`=unlist(subList[["Undetermined"]]), ##n=0
             Inter=unlist(subList[["Inter"]])),
         col=oranges8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", log="", las=2)
dev.off()

##################################################
## for (nbest in 6:10) {
load(sprintf("../external_data/distList_neftel_%d.RData", nbest))
load(sprintf("../external_data/dtList_neftel_%d.RData", nbest))
load(sprintf("../external_data/dtList2_neftel_%d.RData", nbest))


#### Glioma subtypes

## find most divergent cell pair in each subtype
mostDist <- c()
for (disttype in names(dtList2)) {
    both <- unlist(strsplit(disttype, ":"))
    if (both[1] == both[2]) {
        print(both[1])
        l <- unlist(dtList2[[disttype]])
        mostDist[both[1]] = names(l)[which(l==max(l))[1]]
    }
}

represMat <- c()
for (t in names(mostDist)) {
  tidx <- which(subtype == t)
  tmat <- clMat[,tidx]
  ovprof <- apply(tmat, 1, mean)
  tmp <- unlist(strsplit(mostDist[t], ":"))
  represMat <- cbind(represMat, ovprof, clMat[,tmp[2]], clMat[,tmp[3]])
  colnames(represMat)[(ncol(represMat)-2):ncol(represMat)] <- paste(t, c("Overall", "Outlier1", "Outlier2"), sep="_")
}

library(beanplot)
library(beeswarm)
library(viridis)
library(RColorBrewer)
library(gplots)
pdf(sprintf("../plots/neftel_gliosubtype_profiles_%d.pdf", nbest))
##heatmap.2(represMat, col=bluered(100), dendrogram="none", Colv=NA, Rowv=NA, trace="none", mar=c(10,5),
##          scale="none")

span1=1/(ncol(represMat)-1)
span2=1/(nrow(represMat)-1)
par(mar=c(10,4,1,1))
image(t(represMat)[,nrow(represMat):1], col=bluered(100), zlim=c(0,1), axes=F)
axis(1, at=seq(0, 1, span1), label=colnames(represMat), las=2, tick=F, line=NA)
axis(2, at=seq(0, 1, span2), label=paste("CL", nrow(represMat):1, sep=""), las=2, tick=F, line=NA)

image(data.matrix(1:100), col=bluered(100), yaxt="n")
dev.off()
