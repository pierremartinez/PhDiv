tpmsets <- c("filbin", "li", "tirosh_m", "tirosh_o", "venteicher", "neftel")
names(tpmsets) <- c("filbin", "li", "tirosh_m", "tirosh_o", "venteicher", "neftel")

library(multtest)
stars <- function(p) {
  if (p < 0.001)
    return("***")
  if (p < 0.01)
    return("**")
  if (p < 0.05)
    return("*")
  return("")
}

sigstarsmeans <- function(s, m) {
  pv <- c()
  for (set in s) {
    this <- others <- c()
    thissamps <- grep(set, names(m))
    wt <- wilcox.test(m[thissamps], m[-thissamps])
    pv[set] <- wt$p.value
  }
  mt <- mt.rawp2adjp(pv, "Bonferroni")
  corv <- mt$adjp[order(mt$index),"Bonferroni"]
  return(unlist(lapply(corv, stars)))
}

#threshv <- c(2)
threshv <- c(0,1,2,3,5)
#nbestv <- c(8)
nbestv <- c(6:10)
pcaL <- clL <- list()
pdf("../plots/all_samples_all_sets.pdf", width=14)
for (j in 1:length(threshv)) {
  thresh=threshv[j]
  nbest=nbestv[j]
  clList <- pcaList <- list()
  print(c(nbest, thresh))
  for (i in 1:length(tpmsets)) {
    s1 <- tpmsets[i]
    s2 <- names(tpmsets)[i]
    print(s1)
    tmp <- load(sprintf("../external_data/distList_%s_%.1f.RData", s2, thresh))
    pcaList[[s1]] <- get(tmp)
    tmp <- load(sprintf("../external_data/distList_%s_%d.RData", s2, nbest))
    clList[[s1]] <- get(tmp)
  }

  npats <- unlist(lapply(pcaList, function(x){length(x[["Malignant"]])}))
  nobs <- unlist(lapply(pcaList, function(x){unlist(lapply(x[["Malignant"]], length))}))
  samples <- unlist(lapply(pcaList, function(x){names(x[["Malignant"]])}))
  sets <- rep(tpmsets, npats)
  ids <- paste(sets, samples, sep="_")
  names(sets) <- names(ids) <- samples

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
  setcol <- c(purples8[5], greens8[5], blues8[5], reds8[5], greys8[5], oranges8[5])
  names(setcol) <- tpmsets

  cl <- pca <- list()
  for (i in which(nobs > 0)) {
    set=sets[i]
    samp=samples[i]
    id=ids[i]
    cl[[id]] <- clList[[set]][["Malignant"]][[samp]]
    pca[[id]] <- pcaList[[set]][["Malignant"]][[samp]]
  }
  pcaL[[as.character(thresh)]] <- pca
  clL[[as.character(nbest)]] <- cl
  mpca <- unlist(lapply(pca, median))
  mcl <- unlist(lapply(cl, median))

  sigmpca <- sigstarsmeans(s=tpmsets, m=mpca)
  sigmcl <- sigstarsmeans(s=tpmsets, m=mcl)

  par(mar=c(5, 4, 2.5, .5), mfrow=c(2,1))
  boxplot(pca, outline=F, at=order(order(mpca)), las=2, col=setcol[sets[which(nobs>0)]], ylab="Malignant distance",
          main=sprintf("PCA-based (threshold=%d%%)", thresh), cex.axis=0.75, names=samples[which(nobs>0)])
  legend("topleft", fill=setcol, legend=paste(tpmsets, sigmpca), cex=0.75)
  boxplot(cl, outline=F, at=order(order(mcl)), las=2, col=setcol[sets[which(nobs>0)]], ylab="Malignant distance",
          main=sprintf("Cluster-based (n=%d)", nbest), cex.axis=0.75, names=samples[which(nobs>0)])
  legend("topleft", fill=setcol, legend=paste(tpmsets, sigmcl), cex=0.75)
}
dev.off()

library(beeswarm)
pdf("../plots/all_variances.pdf")
for (i in as.character(threshv)) {
  cv <- c()
  pca <- pcaL[[i]]
  mdist <- unlist(lapply(pca, mean))
  setL <- lapply(tpmsets, function(s){mdist[which(sets[nobs>0] == s)]})
  for (set in names(setL)) {
    smean <- mean(setL[[set]])
    ssd <- sd(setL[[set]])
    cv[set] <- ssd / smean
  }
  o <- order(cv)
  o2 <- order(o)
  boxplot(setL, at=o2, outline=F, ylim=c(min(unlist(setL)), max(unlist(setL))), ylab="Mean malignant distance",
          main=sprintf("PCA (threshold=%s%%)", i), names=tpmsets)
  beeswarm(setL, at=o2, #unlist(lapply(setL, length)),
           pwcol=rep(setcol[tpmsets], unlist(lapply(setL, length))), add=T, pch=16)
  boxplot(cv, outline=F, ylim=c(min(cv), max(cv)), ylab="Coefficient of variation", main=sprintf("PCA (threshold=%s%%)", i))
  beeswarm(cv, pwcol=setcol, pch=16, add=T, cex=2)
  legend("topleft", legend=tpmsets, col=setcol, pch=16)
}
for (j in as.character(nbestv)) {
  cv <- c()
  cl <- clL[[j]]
  mdist <- unlist(lapply(cl, mean))
  setL <- lapply(tpmsets, function(s){mdist[which(sets[nobs>0] == s)]})
  for (set in names(setL)) {
    smean <- mean(setL[[set]])
    ssd <- sd(setL[[set]])
    cv[set] <- ssd / smean
  }
  o <- order(cv)
  o2 <- order(o)
  boxplot(setL, at=o2, outline=F, ylim=c(min(unlist(setL)), max(unlist(setL))), ylab="Mean malignant distance",
          main=sprintf("Clusters (n=%s)", j), names=tpmsets)
  beeswarm(setL, at=o2, #unlist(lapply(setL, length)),
           pwcol=rep(setcol[tpmsets], unlist(lapply(setL, length))), add=T, pch=16)
  boxplot(cv, outline=F, ylim=c(min(cv), max(cv)), ylab="Coefficient of variation", main=sprintf("Clusters (n=%s)", j))
  beeswarm(cv, pwcol=setcol, pch=16, add=T, cex=2)
  legend("topleft", legend=tpmsets, col=setcol, pch=16)
}
dev.off()

pcaClMat <- matrix(nc=4, nr=0)
for (i in as.character(threshv)) {
  pca <- pcaL[[i]]
  ##mpca <- unlist(lapply(pca, median))
  print(paste("thresh =",i))
  for (j in as.character(nbestv)) {
    print(paste("nbest =",j))
    cl <- clL[[j]]
    ##mcl <- unlist(lapply(cl, median))
    png(sprintf("../plots/cor_pca_%s_cluster_%s.png", i, j))
    par(mar=c(6.5, 4, 2.5, .5), mfrow=c(1,1))
    ct <- cor.test(unlist(pca), unlist(cl), method="spearman")
    plotcol <- rep(tpmsets, unlist(lapply(tpmsets, function(x){sum(nobs[which(sets==x)])})))
    plot(unlist(pca), unlist(cl), main=sprintf("%s clusters, >%s%% threshold", i, j),
         xlab=paste("PCA-based distance"), ylab="Cluster-based distance", pch=".", cex=2, col=setcol[plotcol], 
       sub=sprintf("tau=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
    abline(lm(unlist(cl) ~ unlist(pca)))
    legend("topleft", col=setcol, legend=tpmsets, pch=16, cex=1)
    pcaL[[as.character(thresh)]] <- pca
    clL[[as.character(nbest)]] <- cl
    dev.off()
    pcaClMat <- rbind(pcaClMat, c(i, j, sprintf("%.2f", ct$estimate), sprintf("%.1e", ct$p.value)))
  }
}

pcaPwMat <- matrix(nc=4, nr=0)
for (i in 1:(length(threshv)-1)) {
  t1=threshv[i]
  pca1 <- pcaL[[as.character(t1)]]
  print(paste("t1 =",t1))
  for (j in (i+1):length(threshv)) {
    t2=threshv[j]
    pca2 <- pcaL[[as.character(t2)]]
    print(paste("  t2 =",t2))
    png(sprintf("../plots/cor_pca_%s_pca_%s.png", t1, t2))
    par(mar=c(6.5, 4, 2.5, .5), mfrow=c(1,1))
    ct <- cor.test(unlist(pca1), unlist(pca2), method="spearman")
    plotcol <- rep(tpmsets, unlist(lapply(tpmsets, function(x){sum(nobs[which(sets==x)])})))
    plot(unlist(pca1), unlist(pca2), main=sprintf("PCA-based distance", i, j),
         xlab=sprintf("Threshold = %d%%", t1), ylab=sprintf("Threshold = %d%%", t2), pch=".", cex=2, col=setcol[plotcol], 
       sub=sprintf("rho=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
    abline(lm(unlist(pca2) ~ unlist(pca1)))
    legend("topleft", col=setcol, legend=tpmsets, pch=16, cex=1)
    dev.off()
    pcaPwMat <- rbind(pcaPwMat, c(t1, t2, sprintf("%.2f", ct$estimate), sprintf("%.1e", ct$p.value)))
  }
}

clPwMat <- matrix(nc=4, nr=0)
for (i in 1:(length(nbestv)-1)) {
  t1=nbestv[i]
  cl1 <- clL[[as.character(t1)]]
  print(paste("t1 =",t1))
  for (j in (i+1):length(nbestv)) {
    t2=nbestv[j]
    cl2 <- clL[[as.character(t2)]]
    print(paste("  t2 =",t2))
    png(sprintf("../plots/cor_cluster_%s_cluster_%s.png", t1, t2))
    par(mar=c(6.5, 4, 2.5, .5), mfrow=c(1,1))
    ct <- cor.test(unlist(cl1), unlist(cl2), method="spearman")
    plotcol <- rep(tpmsets, unlist(lapply(tpmsets, function(x){sum(nobs[which(sets==x)])})))
    plot(unlist(cl1), unlist(cl2), main=sprintf("Cluster-based distance", i, j),
         xlab=sprintf("n = %d", t1), ylab=sprintf("n = %d", t2), pch=".", cex=2, col=setcol[plotcol], 
       sub=sprintf("rho=%.2f; p=%.3f (Spearman correlation)", ct$estimate, ct$p.value))
    abline(lm(unlist(cl2) ~ unlist(cl1)))
    legend("topleft", col=setcol, legend=tpmsets, pch=16, cex=1)
    dev.off()
    clPwMat <- rbind(clPwMat, c(t1, t2, sprintf("%.2f", ct$estimate), sprintf("%.1e", ct$p.value)))
  }
}

colnames(pcaClMat) <- colnames(pcaPwMat) <- colnames(clPwMat) <- c("i1", "i2", "rho", "p")
write.table(pcaClMat, file="../pca_cluster_correlations.txt", quote=F, sep="\t", row.names=F)
write.table(pcaPwMat, file="../pca_pairwise_correlations.txt", quote=F, sep="\t", row.names=F)
write.table(clPwMat, file="../cluster_pairwise_correlations.txt", quote=F, sep="\t", row.names=F)
