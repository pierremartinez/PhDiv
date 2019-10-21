rm(list=ls())
author="tirosh_m"
load("../preprocessed/data_score_AUC_TPM_tirosh_m.Rdata")
load(sprintf("../preprocessed/AUCell_cluster_list_%s.RData", author))

tirosh <- read.table("../external_data/GSE72056_melanoma_single_cell_revised_v2.txt", sep="\t", stringsAsFactors=F, header=T)
meta <- data.matrix(tirosh[1:3,-1])
rownames(meta) <- tirosh[1:3,1]
rownames(meta)[2] <- "malignant" #malignant(1=no,2=yes,0=unresolved) 
rownames(meta)[3] <- "nmType" #(1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)
grp <- meta[3,]
grp[meta[2,]==2] <- 7
grp[meta[2,]==0 | grp==0] <- 8
equivgrp <- c("Tcell", "Bcell", "Macrophage", "Endothelial", "Fibroblast", "NK", "Malignant", "NA")
grp[grp==1]
rm(tirosh)
gc()

cellcat <- meta["malignant",]
allPat <- as.character(sort(unique(meta["tumor",])))

##ind <- get_pca_ind(metapca)

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
distList <- list(m=list(), b=list(), u=list(), inter=list()) ## malign/malign, benign/benign distances
bList <- list() ##intra-benign distances
detailList <- dtList <- list()
for (p in allPat) {
  print(p)
  pmat <- clMat[,which(meta["tumor",] == p)]
  combs <- combn(colnames(pmat), 2)
  pdist <- dist(t(pmat), "euclidean")
  distcat <- unlist(mclapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    if (meta["malignant", c1] != meta["malignant", c2])
      return("inter")
    if (meta["malignant", c1] == 0) ##unresolved
      return("u")
    if (meta["malignant", c1] == 1) ##benign
      return("b")
    return("m")
  }, mc.cores=8))
  distname <- unlist(mclapply(1:ncol(combs), function(i){
      paste(c(p,combs[,i]), collapse=":")
  }, mc.cores=8))
  names(pdist) <- distname
  distList[["Malignant"]][[p]] <- pdist[which(distcat == "m")]
  distList[["Benign"]][[p]] <- pdist[which(distcat == "b")]
  distList[["Undefined"]][[p]] <- pdist[which(distcat == "u")]
  distList[["Inter"]][[p]] <- pdist[which(distcat == "inter")]

  ## benign only
  bidx <- which(distcat=="b")
  bcat <- unlist(mclapply(bidx, function(i){ ##(1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK
    c1=combs[1,i]
    c2=combs[2,i]
    if (meta["nmType", c1] != meta["nmType", c2])
      return("inter")
    if (meta["nmType", c1] == 1)
      return("T")
    if (meta["nmType", c1] == 2)
      return("B")
    if (meta["nmType", c1] == 3)
      return("Macro")
    if (meta["nmType", c1] == 4)
      return("Endo")
    if (meta["nmType", c1] == 5)
      return("CAF")
    return("NK")
  }, mc.cores=8))
  bdist <- pdist[bidx]
  names(bdist) <- bcat
  bList[["Tcell"]][[p]] <- bdist[which(bcat == "T")]
  bList[["Bcell"]][[p]] <- bdist[which(bcat == "B")]
  bList[["Macrophage"]][[p]] <- bdist[which(bcat == "Macro")]
  bList[["Fibroblast"]][[p]] <- bdist[which(bcat == "CAF")]  
  bList[["Endothelial"]][[p]] <- bdist[which(bcat == "Endo")]  
  bList[["NK"]][[p]] <- bdist[which(bcat == "NK")]  
  bList[["inter"]][[p]] <- bdist[which(bcat == "inter")]  

  detailcat <- unlist(mclapply(1:ncol(combs), function(i){
    c1=combs[1,i]
    c2=combs[2,i]
    return(paste(sort(c(equivgrp[grp[c1]], equivgrp[grp[c2]])), collapse=":"))
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
grpo <- unlist(lapply(1:length(unique(grp)), function(x){which(grp==x)}))

save(distList, file=sprintf("../external_data/distList_tirosh_m_%d.RData", nbest))
save(bList, file=sprintf("../external_data/bList_tirosh_m_%d.RData", nbest))
save(dtList, file=sprintf("../external_data/dtList_tirosh_m_%d.RData", nbest))

## for (nbest in 6:10) {
load(sprintf("../external_data/distList_tirosh_m_%d.RData", nbest))
load(sprintf("../external_data/bList_tirosh_m_%d.RData", nbest))
load(sprintf("../external_data/dtList_tirosh_m_%d.RData", nbest))

meanm <- unlist(lapply(distList[["Malignant"]], mean, na.rm=T))
meanb <- unlist(lapply(distList[["Benign"]], mean, na.rm=T))
library(vegan)
simp <- unlist(lapply(allPat, function(p){
  diversity(table(meta["nmType",which(meta["tumor",] == p & meta["malignant",] == 1)]), index="simpson")
}))
shan <- unlist(lapply(allPat, function(p){
  diversity(table(meta["nmType",which(meta["tumor",] == p & meta["malignant",] == 1)]))
}))


cellorder <- c("Tcell", "Bcell", "Macrophage", "Fibroblast", "Endothelial", "NK", "Malignant", "NA")
patgrpmat <- rbind(equivgrp[grp], meta["tumor",])
patgrp <- 1:length(c(cellorder, allPat))
names(patgrp) <- c(cellorder, allPat)
pgm2 <- matrix(NA, nc=ncol(patgrpmat), nr=nrow(patgrpmat))
pgm2[1,] <- patgrp[patgrpmat[1,]]
pgm2[2,] <- patgrp[patgrpmat[2,]]
grpcol <- c(rainbow, viridis(length(allPat)))
names(grpcol) <- names(patgrp)
patgrpo <- order(pgm2[1,], partial=pgm2[2,])

pdf(sprintf("../plots/tirosh_melanoma_distances_%d.pdf", nbest))
heatmap.2(clMat[,patgrpo], Colv=F, Rowv=T, col=bluered(100), trace="none", labCol=NA, dendrogram="row",
          labRow=paste("CL", 1:nbest, sep=""), cexRow=.6, mar=c(1,3))
par(mar=c(.5,.5,.5,.5), mfrow=c(2,1))
image(t(pgm2[2:1,patgrpo]), col=grpcol, axes=F)
plot.new()
legend("topleft", legend=cellorder, fill=grpcol[1:length(cellorder)], horiz=F, cex=1, box.lty=0)
legend("topright", legend=allPat, fill=grpcol[(1+length(cellorder):length(grpcol))], horiz=F, cex=.75, box.lty=0)

ct <- cor.test(simp, meanb, method="spearman")
plot(simp, meanb, ylab="Mean malignant-malignant distance", xlab="Simpson diversity (Benign compartment)", col=blues8[5], pch=16,
     sub=sprintf("tau=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
abline(lm(meanb ~ simp))
ct <- cor.test(shan, meanb, method="spearman")
plot(shan, meanb, ylab="Mean malignant-malignant distance", xlab="Shannon diversity (Benign compartment)", col=blues8[5], pch=16,
     sub=sprintf("tau=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
abline(lm(meanb ~ shan))

beanplot(list(Malignant=unlist(distList[["Malignant"]]),
              Tcell=unlist(bList[["Tcell"]]),
              Bcell=unlist(bList[["Bcell"]]),
              Macrophage=unlist(bList[["Macrophage"]]),
              Fibroblast=unlist(bList[["Fibroblast"]]),
              Endothelial=unlist(bList[["Endothelial"]]),
              NK=unlist(bList[["NK"]]),
              Undefined=unlist(distList[["Undefined"]]),
              Inter=unlist(distList[["Inter"]])),
         col=blues8[c(5,1,3,8)], cut=.1, what=c(1,1,1,0), ylab="Distance", las=2, log="")
dev.off()

## for (nbest in 6:10) {
load(sprintf("../external_data/distList_tirosh_m_%d.RData", nbest))
load(sprintf("../external_data/bList_tirosh_m_%d.RData", nbest))
load(sprintf("../external_data/dtList_tirosh_m_%d.RData", nbest))

## find most divergent cell pair in each type
mostDist <- c()
for (disttype in names(dtList)) {
    both <- unlist(strsplit(disttype, ":"))
    if (both[1] == both[2]) {
        print(both[1])
        l <- unlist(dtList[[disttype]])
        mostDist[both[1]] = names(l)[which(l==max(l))[1]]
    }
}

represMat <- c()
for (t in names(mostDist)) {
    tidx <- which(equivgrp[grp] == t)
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
pdf(sprintf("../plots/tirosh_melanoma_subtype_profiles_%d.pdf", nbest))
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
