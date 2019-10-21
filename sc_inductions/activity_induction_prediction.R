rm(list=ls())

library(gdata) ; library(gplots) ; library(BPSC) ; library(multtest) ; library(beeswarm) ; library(ggplot2)
library(stringr) ; library(pROC) ; library(gdata) 

############################################################################################################
############################################################################################################
############################################################################################################


################# Importation des données

genedat <- read.table("emt_dnarep_glyco_induction_sc.csv", header=T, fill=T, stringsAsFactors=F, sep="\t")
c1dat <- read.xls("Biomark_output_table.xls", stringsAsFactors=F, header=F, sheet=3)

colnames(c1dat) <- c("Well", "Cell", "Type", "rConc", "Gene", "Type", "Value", "Quality", "Call", "Threshold", "Tm_In", "Tm_Out", "Peak.Ratio")
c1dat <- c1dat[1:(48*48),]
c1dat$Gene[c1dat$Gene == "SCL16A3"] <- "SLC16A3"
c1dat$Gene[c1dat$Gene == "CTFG"] <- "CTGF"
c1dat$Gene[c1dat$Gene == "R.spk1mid"] <- "spike1"
c1dat$Gene[c1dat$Gene == "R.spk4mid"] <- "spike4"
genes <- intersect(genedat$Genes, unique(c1dat$Gene))

c1dat$Cell <- gsub("^1_1", "DNA_repair", c1dat$Cell)
c1dat$Cell <- gsub("^1_2", "Control", c1dat$Cell)
c1dat$Cell <- gsub("^1_3", "Glyco", c1dat$Cell)
c1dat$Cell <- gsub("^1_4", "EMT", c1dat$Cell)
cells <- unique(c1dat$Cell)
cells <- c(cells[grep("EMT", cells)], cells[grep("DNA_repair", cells)], cells[grep("Glyco", cells)], cells[grep("Control", cells)], "C-", "C+")

### 46 gènes (+2 spikes) et 48 cellules


################# Filtre et matrice Ct

failMat <- ctMat <- matrix(NA, nr=length(genes)+2, nc=length(cells))
colnames(ctMat) <- colnames(failMat) <- cells
rownames(ctMat) <- rownames(failMat) <- c(genes, "spike1", "spike4")
for (i in 1:nrow(c1dat)) {
  ctMat[c1dat$Gene[i],c1dat$Cell[i]] <- as.numeric(c1dat$Value[i])
  failMat[c1dat$Gene[i],c1dat$Cell[i]] <- ifelse(c1dat$Call[i] == "Pass", 0, 1)
}
cellfail <- apply(failMat, 2, sum)/nrow(failMat)
genefail <- apply(failMat, 1, sum)/ncol(failMat)

cellthresh=0.3
genethresh=0.3

ctMatThresh <- ctMat
ctMatThresh[which(genefail > genethresh),] <- 999
ctMatThresh[,which(cellfail > cellthresh)] <- 999

##ctMat : matrice 48x48 , en ligne les gènes + 2 spikes et en colonne les cellules
##ctMatThresh : matrice 48x48 , en ligne les gènes + 2 spikes et en colonne les cellules, plus filtre



################# Remplacer 999 par NA

remplir_NA=function(m){
  m1=m
  for(i in 1:nrow(m)){
    ind=which(m[i,]==999)
    m1[i,ind]=NA
  }
  return(m1)
}
ctMatThresh_NA=remplir_NA(ctMatThresh)


data.Ct_tresh=ctMatThresh_NA[-c(which(genefail > genethresh)),-c(which(cellfail > genethresh))]
### Matrice Ct sans les gènes et cellules qui ont fail
### Matrice 36x37



############################################################################################################
############################       Etape de pré-normalisation       ####################################
############################################################################################################


### Récupération des spikes
spike1.thresh=data.Ct_tresh["spike1",]
spike4.thresh=data.Ct_tresh["spike4",]

### Test de la normalité
shapiro.test(spike1.thresh) ; shapiro.test(spike4.thresh)


### Enlever les valeurs extrêmes
remove.cell.spike=function(v_spike){
  while(shapiro.test(v_spike)$p.value<0.05){
    dist=abs(v_spike-summary(v_spike)[3])
    ind=which(dist==max(dist,na.rm=T))
    v_spike[ind]=NA
  }
  return(v_spike)
}

spike1.thresh.r=remove.cell.spike(spike1.thresh)
spike4.thresh.r=remove.cell.spike(spike4.thresh)


### Test de Shapiro pour vérifier la normalité
shapiro.test(spike1.thresh.r) ; shapiro.test(spike4.thresh.r)

which(is.na(spike1.thresh.r)) ; which(is.na(spike4.thresh.r))

spike.thresh.r=cbind(spike1.thresh.r,spike4.thresh.r)
colnames(spike.thresh.r)=c("spike1","spike4")

cell.remove=unique(c(which(is.na(spike1.thresh.r)),which(is.na(spike4.thresh.r))))

data.Ct_tresh.r=data.Ct_tresh[,-c(cell.remove)]
###Données 36x36


###  Moyenne et ecart type
spike.tab.thresh=matrix(c(mean(spike1.thresh,na.rm=T),sd(spike1.thresh,na.rm=T),
                          mean(spike4.thresh,na.rm=T),sd(spike4.thresh,na.rm=T)),ncol=2,byrow=T)
colnames(spike.tab.thresh)=c("Moyenne","Ecart-type")
rownames(spike.tab.thresh)=c("spike1","spike4")

spike.tab.thresh.r=matrix(c(mean(spike1.thresh.r,na.rm=T),sd(spike1.thresh.r,na.rm=T),
                            mean(spike4.thresh.r,na.rm=T),sd(spike4.thresh.r,na.rm=T)),ncol=2,byrow=T)
colnames(spike.tab.thresh.r)=c("Moyenne","Ecart-type")
rownames(spike.tab.thresh.r)=c("spike1","spike4")

###  Liste avec les tableaux avant/après
spike.mean.sd.thresh=list("Avant"=spike.tab.thresh,"Après"=spike.tab.thresh.r)



############################################################################################################
###################################        Normalisation       ####################################
############################################################################################################


ajustement=function(spike){
  m=mean(spike,na.rm=T)
  v=NULL
  for(i in 1:length(spike)){
    v[i]=m-spike[i]
  }
  names(v)=names(spike)
  return(v)
}

adj.spike1.thresh=ajustement(spike1.thresh.r)
adj.spike4.thresh=ajustement(spike4.thresh.r)

normalisation=function(m,adj1,adj2){
  m1=matrix(ncol=ncol(m),nrow=nrow(m))
  colnames(m1)=colnames(m)
  rownames(m1)=rownames(m)
  for(j in 1:ncol(m)){
    m1[,j]=m[,j]+mean(c(adj1[j],adj2[j]),na.rm=T)
  }
  return(m1)
}


### Données normalisées
data.Ct_tresh.r.norm=normalisation(data.Ct_tresh.r,adj.spike1.thresh,adj.spike4.thresh)


### Spikes normalisés
spike1.thresh.r.norm=data.Ct_tresh.r.norm["spike1",]
spike4.thresh.r.norm=data.Ct_tresh.r.norm["spike4",]


spike.tab.thresh.r.norm=matrix(c(mean(spike1.thresh.r.norm,na.rm=T),sd(spike1.thresh.r.norm,na.rm=T),
                                 mean(spike4.thresh.r.norm,na.rm=T),sd(spike4.thresh.r.norm,na.rm=T)),ncol=2,byrow=T)
colnames(spike.tab.thresh.r.norm)=c("Moyenne","Ecart-type")
rownames(spike.tab.thresh.r.norm)=c("spike1","spike4")


### Moyenne et écart-type après correction
spike.tab.thresh.r.norm


############################################################################################################
###################################        Inférence        #######################################
############################################################################################################

##### Transformation Cq en RNAm
matrice_RNA=function(m1){
  m=matrix(ncol=ncol(m1),nrow=nrow(m1))
  rownames(m)=rownames(m1)
  colnames(m)=colnames(m1)
  for(j in 1:ncol(m1)){
    for (i in 1:nrow(m1)){
      
      m[i,j]=ifelse(is.na(m1[i,j])==TRUE, 0, 48 * 45 * 2 ^ (30 - 18 - m1[i,j]))
    }
  }
  return(m)
}


## Matrice avec les valeurs de transcrits
data.thresh=matrice_RNA(data.Ct_tresh.r.norm)


##### Gene à supprimer

gene.EMT=genedat$Genes[grep("EMT",genedat$GeneSet)]
gene.DNA_repair=genedat$Genes[grep("DNA_repair",genedat$GeneSet)]
gene.glyco=genedat$Genes[grep("Glyco",genedat$GeneSet)]


names(which(genefail > genethresh))

sapply(names(which(genefail > genethresh)),function(x){ grep(x,gene.glyco)})
sapply(names(which(genefail > genethresh)),function(x){ grep(x,gene.EMT)})
sapply(names(which(genefail > genethresh)),function(x){ grep(x,gene.DNA_repair)})

gene.EMT.thresh=gene.EMT[-c(3,4,6,10,11,13)]
gene.DNA_repair.thresh=gene.DNA_repair[-c(1,2,7,8,9,10)]


##### Données par activités
data.glyco.thresh=data.thresh[gene.glyco,]
data.EMT.thresh=data.thresh[gene.EMT.thresh,]
data.DNA_repair.thresh=data.thresh[gene.DNA_repair.thresh,]

##Matrice de taille i x 36, où i est le nombre de gènes par activité

#####################################################################
####            Pierre addition: heatmap per activity            ####
#####################################################################
library(gplots)
library(BPSC)
library(multtest)
condequiv <- c("EMT", "DNA_repair", "Glycolysis")
names(condequiv) <- c("EMT", "DNA_repair", "Glyco")
pdf("expression_activ_VS_normal.pdf")
for (cond in names(condequiv)) {
  ce <- condequiv[cond]
  ceg <- intersect(rownames(data.Ct_tresh.r.norm), genedat$Genes[grep(ce, genedat$GeneSet)])
  cec <- colnames(data.Ct_tresh.r.norm)[c(grep(cond, colnames(data.Ct_tresh.r.norm)), grep("Control", colnames(data.Ct_tresh.r.norm)))]
  ceco <- setdiff(colnames(data.Ct_tresh.r.norm), c(cec, "C+"))
  delim <- grep("Control", cec)[1]-1
  ct <- data.Ct_tresh.r[ceg,cec]
  ## if (cond == "EMT") {
  ##   ceg=c(ceg,"ZEB1")
  ## }
  ##ct[ct==999] <- NA
  nct <- data.Ct_tresh.r.norm[ceg,cec]
  ##nct[is.na(ct)] <- NA
  norm <- data.thresh[ceg,cec]
  norm[is.na(ct)] <- NA
  bignorm <- data.thresh[ceg,c(cec,ceco)]
  bignorm[is.na(data.Ct_tresh.r[ceg,c(cec,ceco)])] <- NA
  heatmap.2(ct, Rowv=F, Colv=F, dendrogram="none", trace="none", scale="row", col=bluered(25), margin=c(8,8),
            main=paste(ce, "Raw Ct"), colsep=delim, sepcol="black", na.color="black")
  heatmap.2(nct, Rowv=F, Colv=F, dendrogram="none", trace="none", scale="row", col=bluered(25), margin=c(8,8),
            main=paste(ce, "Normalised Ct"), colsep=delim, sepcol="black", na.color="black")
  heatmap.2(norm, Rowv=F, Colv=F, dendrogram="none", trace="none", scale="row", col=bluered(25), margin=c(8,8),
            main=paste(ce, "Normalised RNA molecules"), colsep=delim, sepcol="black", na.color="black")
  heatmap.2(bignorm, Rowv=F, Colv=F, dendrogram="none", trace="none", scale="row", col=bluered(25), margin=c(8,8),
            main=paste(ce, "Normalised RNA molecules"), colsep=delim, sepcol="black", na.color="black")
}
dev.off()

###############################   Distribution Beta-Poisson   #######################################

### Modélisation Beta-Poisson sur une ligne
estim_BP=function(x){
  t=estimateBP(x[which(!is.na(x))], para.num=4, tbreak.num=20)
  return(t)
}

testFunction <- function (x) {
  return(tryCatch(estim_BP(x), error=function(e) NA))
}

### Modélisation Beta-Poisson sur une matrice
BP_function=function(m,cond){
  ind=grep(cond,colnames(m))
  l=list()
  for(i in 1:nrow(m)){
    l[[rownames(m)[i]]]=list(testFunction(m[i,ind]),testFunction(m[i,-ind]))
    names(l[[rownames(m)[i]]])=c(cond,"Control")
  }
  return(l)
}

BP_glyco.thresh=BP_function(data.glyco.thresh,"Glyco")
BP_EMT.thresh=BP_function(data.EMT.thresh,"EMT")
BP_DNA_repair.thresh=BP_function(data.DNA_repair.thresh,"DNA_repair")


### Récupération des coeffs
coeff.BP=function(l,cond){
  m1=matrix(ncol=4,nrow=length(l))
  m2=matrix(ncol=4,nrow=length(l))
  rownames(m1)=rownames(m2)=names(l)
  colnames(m1)=colnames(m2)=c("alpha","beta","lambda1","lambda2")
  for(i in 1:length(l)){
    if(class(l[[i]][[1]])=="list"){
      m1[i,]=c(l[[i]][[1]]$par)
      m2[i,]=c(l[[i]][[2]]$par)
    }else{
      m1[i,]=m2[i,]=rep(NA,4)
    }
  }
  l2=list(m1,m2)
  names(l2)=c(cond,"Control")
  return(l2)
}

coeff_BP_glyco.thresh=coeff.BP(BP_glyco.thresh,"Glycolysis")
coeff_BP_EMT.thresh=coeff.BP(BP_EMT.thresh,"EMT")
coeff_BP_DNA_repair.thresh=coeff.BP(BP_DNA_repair.thresh,"DNA_repair")



#############################    Intervalle de confiance   #######################################

### Ecart-type pour construire les intervalles
sd.spike.thresh=mean(spike.tab.thresh.r.norm[,2])

## Intervalle pour une valeur de Ct
interval_Ct.2=function(valCt,sd){
  IC=cbind(valCt,valCt-sd,valCt+sd)
  colnames(IC)=c("Valeur observée","Borne inf","Borne sup")
  return(IC)
}

## Sans la première colonne
interval_Ct=function(valCt,sd){
  IC=cbind(valCt-sd,valCt+sd)
  colnames(IC)=c("Borne inf","Borne sup")
  return(IC)
}

## Ct transformée en nombre de transcrit
Ct_RNA=function(val){
  RNAm=48*45*2^(30-18-val)
  return(RNAm)
}

## Nombre de transcrit transformé en Ct
RNA_Ct=function(x){
  Ct=(-log(x/(48*45))/log(2)) + 30 - 18
  return(Ct)
}

### Intervalle pour une valeur de RNA
intervals_vRNA=function(v,sd){
  mat=matrix(ncol=3,nrow=length(v))
  mat[,1]=v
  Ct_v=RNA_Ct(v)
  mat[,c(3,2)]=Ct_RNA(interval_Ct(Ct_v,sd))
  colnames(mat)=c("Valeur","Borne inf","Borne sup")
  return(mat)
}

intervals=function(m,cond,sd,gene){
  ind=grep(cond,colnames(m))
  m1=m[,ind] ; m2=m[,-c(ind)]
  
  IC1=matrix(nrow=length(m1[gene,]),ncol=5)
  IC2=matrix(nrow=length(m2[gene,]),ncol=5)
  
  IC1=cbind(m1[gene,],intervals_vRNA(m1[gene,],sd)[,c(2,3)])
  IC2=cbind(m2[gene,],intervals_vRNA(m2[gene,],sd)[,c(2,3)])
  colnames(IC1)=colnames(IC2)=c("RNAm","Borne inf","Borne sup")
  l=list(IC1,IC2)
  names(l)=c(cond,"Control")
  return(l)
}


### Tous les intervalles
intervals.all=function(m,cond,sd,listgene){
  l=list()
  for(i in 1:length(listgene)){
    l[[listgene[i]]]=intervals(m,cond,sd,listgene[i])
  }
  return(l)
}

IC_glyco.thresh=intervals.all(data.glyco.thresh,"Glyco",sd.spike.thresh,gene.glyco)
IC_EMT.thresh=intervals.all(data.EMT.thresh,"EMT",sd.spike.thresh,gene.EMT.thresh)
IC_DNA_repair.thresh=intervals.all(data.DNA_repair.thresh,"DNA_repair",sd.spike.thresh,gene.DNA_repair.thresh)


###On change les intervalles (+-1)

correction.IC=function(IC){
  l=list()
  n1=names(IC[[1]])[1]
  n2=names(IC[[1]])[2]
  
  for(i in 1:length(IC)){
    l[[i]]=list()
    l[[i]][[n1]]=IC[[i]][[1]]
    l[[i]][[n2]]=IC[[i]][[2]]
    
    dist1=l[[i]][[n1]][,3]-l[[i]][[n1]][,2]
    dist2=l[[i]][[n2]][,3]-l[[i]][[n2]][,2]
    
    ind1=which(dist1<1)
    ind2=which(dist2<1)
    
    l[[i]][[n1]][ind1,c(2,3)]=c(l[[i]][[n1]][ind1,1]-1,l[[i]][[n1]][ind1,1]+1)
    l[[i]][[n2]][ind2,c(2,3)]=c(l[[i]][[n2]][ind2,1]-1,l[[i]][[n2]][ind2,1]+1)
    
    ind3=which(l[[i]][[n1]][,2]<0) ; ind4=which(l[[i]][[n2]][,2]<0)
    l[[i]][[n1]][ind3,2]=0 ; l[[i]][[n2]][ind4,2]=0
  }
  
  names(l)=names(IC)
  return(l)
}

IC_glyco.thresh.2=correction.IC(IC_glyco.thresh)
IC_EMT.thresh.2=correction.IC(IC_EMT.thresh)
IC_DNA_repair.thresh.2=correction.IC(IC_DNA_repair.thresh)



###Creation d'une matrice enlevant chaque valeur observée
remove_val_obs=function(m){
  mat=matrix(ncol=length(m),nrow=length(m))
  for(i in 1:nrow(mat)){
    mat[i,]=m
    if(is.na(mat[i,i])==TRUE){
      mat[i,]=rep(NA,ncol(mat))
    }else{
      mat[i,i]=NA
    }
  }
  return(mat)
}

m_sans_obs=function(m,cond){
  ind=grep(cond,colnames(m))
  m1=m[,ind]
  m2=m[,-ind]
  
  l1=list()
  for(i in 1:nrow(m1)){
    l1[[i]]=remove_val_obs(m1[i,])
  }
  
  l2=list()
  for(i in 1:nrow(m2)){
    l2[[i]]=remove_val_obs(m2[i,])
  }
  
  names(l1)=names(l2)=rownames(m1)
  
  l=list(l1,l2)
  names(l)=c(cond,"Control")
  
  return(l)
}

sansobs_Glyco.thresh=m_sans_obs(data.glyco.thresh,"Glyco")
sansobs_EMT.thresh=m_sans_obs(data.EMT.thresh,"EMT")
sansobs_DNA_repair.thresh=m_sans_obs(data.DNA_repair.thresh,"DNA_repair")


### Beta-Poisson sans la valeur observée
recup.coeff=function(l){
  m=NULL
  if(class(l)=="list"){
    m=c(l$par)
  }else{
    m=rep(NA,4)
  }
  names(m)=c("alpha","beta","lambda1","lambda2")
  return(m)
}

BP_function_sansobs=function(l){
  
  l1=list() ; l2=list()
  for(i in 1:length(l[[1]])){
    print(names(l[[1]])[i])
    l1[[i]]=list() ; l2[[i]]=list()
    l1[[i]]=apply(l[[1]][[i]],1,testFunction)
    l2[[i]]=apply(l[[2]][[i]],1,testFunction)
  }
  
  names(l1)=names(l2)=names(l[[1]])
  
  coeff1=list() ;   coeff2=list()
  for(j in 1:length(l1)){
    coeff1[[j]]=t(sapply(l1[[j]],recup.coeff))
    coeff2[[j]]=t(sapply(l2[[j]],recup.coeff))
  }
  
  names(coeff1)=names(coeff2)=names(l[[1]])
  coeff=list(coeff1,coeff2)
  names(coeff)=names(l)
  return(coeff)
}

coeff_sansobs_glyco.thresh=BP_function_sansobs(sansobs_Glyco.thresh)
coeff_sansobs_EMT.thresh=BP_function_sansobs(sansobs_EMT.thresh)
coeff_sansobs_DNA_repair.thresh=BP_function_sansobs(sansobs_DNA_repair.thresh)


### Probabilité
proba_rBP=function(IC,coeff,coeff_sansobs){
  proba=list()
  cond=names(coeff)[1]
  cond2=names(coeff)[2]
  for(i in 1:length(IC)){
    if(is.na(coeff[[1]][i,1])==FALSE){
      proba1=matrix(nrow=nrow(IC[[i]][[1]]),ncol=5)
      proba2=matrix(nrow=nrow(IC[[i]][[2]]),ncol=5)
      colnames(proba1)=colnames(proba2)=c("RNA","Borne inf", "Borne sup",paste("Probabilité",cond,sep=" "),
                                          paste("Probabilité",cond2,sep=" "))
      rownames(proba1)=rownames(IC[[i]][[1]])
      rownames(proba2)=rownames(IC[[i]][[2]])
      
      proba1[,c(1,2,3)]=IC[[i]][[1]][,c(1,2,3)]
      proba2[,c(1,2,3)]=IC[[i]][[2]][,c(1,2,3)]
      
      rBP_cond=rBP(10000,coeff[[1]][i,])
      rBP_cond2=rBP(10000,coeff[[2]][i,])
      
      for(j in 1:nrow(proba1)){
        bons=which(rBP_cond2 >= proba1[j,2] & rBP_cond2 <=proba1[j,3])
        x1=rBP(10000,coeff_sansobs[[1]][[i]][j,])
        bons1=which(x1 >= proba1[j,2] & x1 <=proba1[j,3])
        proba1[j,c(4,5)]=c(length(bons1)/length(x1),length(bons)/length(rBP_cond2))
      }
      
      for(k in 1:nrow(proba2)){
        bons=which(rBP_cond >= proba2[k,2] & rBP_cond <=proba2[k,3])
        x2=rBP(10000,coeff_sansobs[[2]][[i]][k,])
        bons2=which(x2 >= proba2[k,2] & x1 <=proba2[k,3])
        proba2[k,c(4,5)]=c(length(bons)/length(rBP_cond),length(bons2)/length(x2))
      }
      
      proba[[i]]=rbind(proba1,proba2)
    }
  }
  names(proba)=names(IC)
  return(proba)
}

proba_glyco.thresh=proba_rBP(IC_glyco.thresh.2,coeff_BP_glyco.thresh,coeff_sansobs_glyco.thresh)
proba_EMT.thresh=proba_rBP(IC_EMT.thresh.2,coeff_BP_EMT.thresh,coeff_sansobs_EMT.thresh)
proba_DNA_repair.thresh=proba_rBP(IC_DNA_repair.thresh.2,coeff_BP_DNA_repair.thresh,coeff_sansobs_DNA_repair.thresh)


### Vraisemblance
likelihood=function(proba){
  l=list()
  cond=substr(colnames(proba[[1]])[4],13,nchar(colnames(proba[[1]])[4]))
  cond2=substr(colnames(proba[[1]])[5],13,nchar(colnames(proba[[1]])[5]))
  
  for(i in 1:length(proba)){
    if(class(proba[[i]])=="matrix"){
      l[[i]]=cbind(proba[[i]],proba[[i]][,4]/(proba[[i]][,4]+proba[[i]][,5]),
                   proba[[i]][,5]/(proba[[i]][,4]+proba[[i]][,5]))
      colnames(l[[i]])=c(colnames(proba[[i]]),paste("Likelihood",cond),paste("Likelihood",cond2))
      ## hack added - Pierre
      l[[i]][,6][is.na(l[[i]][,6]) & l[[i]][,4] == 0] <- 0.5
      l[[i]][,7][is.na(l[[i]][,7]) & l[[i]][,4] == 0] <- 0.5
    }else{
      l[[i]]=NULL
    }
  }
  names(l)=names(proba)
  return(l)
}

likelihood_glyco.thresh=likelihood(proba_glyco.thresh)
likelihood_EMT.thresh=likelihood(proba_EMT.thresh)
likelihood_DNA_repair.thresh=likelihood(proba_DNA_repair.thresh)


### Courbes ROC
calcul.roc.2=function(l,cond){
  roc1=list() ; roc2=list() ; r=list()
  ind=grep("Likelihood",colnames(l[[1]]))
  ind1=grep(cond,rownames(l[[1]]))
  
  n1=length(ind1)
  n2=nrow(l[[1]])-n1
  
  pdf(paste("ROC_curve_",cond,".pdf",sep=""))
  
  for(i in names(l)){
    roc1[[i]]=roc(response= c(rep(1,n1),rep(0,n2)), predictor=l[[i]][,ind[1]])
    roc2[[i]]=roc(response= c(rep(0,n1),rep(1,n2)), predictor=l[[i]][,ind[2]])
    
    par(mar=c(7,5,3,3),xpd=F)
    plot(1-roc1[[i]]$specificities, roc1[[i]]$sensitivities, type="l", 
         xlab="1 - Specificity", ylab="Sensitivity", col="blue",
         main=paste(cond,"from Control",sep=" "),font.sub=4,cex.sub=0.8,cex.axis=0.8)
    points(1-roc2[[i]]$specificities, roc2[[i]]$sensitivities, type="l", col="red")
    points(x=c(0,1), y=c(0,1), type="l", lty=2)
    legend("bottomright",inset=0.05, col=c("blue", "red"), legend=c(cond, "Control"), 
           lty=1,box.lty = 0,horiz=T,cex=0.9)
    legend("topleft",inset=0.07,legend=i,text.font=4,box.lty=0,cex=0.8)
    par(xpd=T)
    text(0.5,-0.25,paste("AUC",cond,"=", round(roc1[[i]]$auc,digits=4),"   ", 
                         "AUC Control =",round(roc2[[i]]$auc,digits=4),sep=" "),
         font=2, cex=0.8)
    
    r[[i]]=list(roc1[[i]],roc2[[i]])
    names(r[[i]])=c(cond,"Control")
  }
  dev.off()
  return(r)
}

roc_glyco.thresh=calcul.roc.2(likelihood_glyco.thresh,"Glyco")
roc_DNA_repair.thresh=calcul.roc.2(likelihood_DNA_repair.thresh,"DNA_repair")
roc_EMT.thresh=calcul.roc.2(likelihood_EMT.thresh,"EMT")


### Récupération des AUC
calcul.auc=function(roc){
  m.auc=matrix(ncol=2,nrow=length(roc))
  rownames(m.auc)=names(roc)
  colnames(m.auc)=names(roc[[1]])
  for(i in 1:length(roc)){
    m.auc[i,]=c(roc[[i]][[1]]$auc,roc[[i]][[2]]$auc)
  }
  return(m.auc)
}

auc_glyco.thresh=calcul.auc(roc_glyco.thresh)
auc_DNA_repair.thresh=calcul.auc(roc_DNA_repair.thresh)
auc_EMT.thresh=calcul.auc(roc_EMT.thresh)

## Combinaison de max 10 gènes
combinaison=function(liste,nbmax){
  l=list()
  for(i in 2:min(nbmax,length(liste))){
    l[[c(paste("Combinaison de",i,"gènes"))]]=combn(liste,i)
  }
  return(l)
}

comb_glyco=combinaison(gene.glyco,10)
comb_EMT=combinaison(gene.EMT.thresh,10)
comb_DNA_repair=combinaison(gene.DNA_repair.thresh,10)


## Noms des combinaisons
m.nom=function(comb){
  v=list()
  for(j in 1:length(comb)){
    v[[j]]=matrix(ncol=1,nrow=ncol(comb[[j]]))
    for(i in 1:ncol(comb[[j]])){
      x=comb[[j]][1,i]
      for(k in 2:nrow(comb[[j]])){
        x=str_c(x,comb[[j]][k,i],sep=", ")
      }
      v[[j]][i,]=x
    }
  }
  names(v)=names(comb)
  return(v)
}

n.comb_glyco=m.nom(comb_glyco)
n.comb_EMT=m.nom(comb_EMT)
n.comb_DNA_repair=m.nom(comb_DNA_repair)

## Liste organisée pour faire l'inférence
list.mcomb=function(proba,comb,ncomb,cond){
  l=list() ; l1=list() ; l2=list()
  
  for(i in 1:length(comb)){
    k=nrow(comb[i])
    l1[[i]]=list() ; l2[[i]]=list() 
    
    for(j in 1:ncol(comb[[i]])){
      n.gene=comb[[i]][,j]
      l1[[i]][[j]]=as.matrix(c(proba[[n.gene[1]]][,6]))
      l2[[i]][[j]]=as.matrix(c(proba[[n.gene[1]]][,7]))
      
      for(z in 2:length(n.gene)){
        l1[[i]][[j]]=cbind(l1[[i]][[j]],c(proba[[n.gene[z]]][,6]))
        l2[[i]][[j]]=cbind(l2[[i]][[j]],c(proba[[n.gene[z]]][,7]))
      }
      colnames(l1[[i]][[j]])=c(n.gene) ; colnames(l2[[i]][[j]])=c(n.gene)
      
    }
    names(l1[[i]])=ncomb[[i]] ; names(l2[[i]])=ncomb[[i]]
  }
  l=list(l1,l2)
  names(l)=c(cond,"Control")
  names(l[[1]])=names(comb) ; names(l[[2]])=names(comb)
  return(l)
}

l.mcomb_glyco=list.mcomb(likelihood_glyco.thresh,comb_glyco,n.comb_glyco,"Glycolysis")
l.mcomb_EMT=list.mcomb(likelihood_EMT.thresh,comb_EMT,n.comb_EMT,"EMT")
l.mcomb_DNA_repair=list.mcomb(likelihood_DNA_repair.thresh,comb_DNA_repair,n.comb_DNA_repair,"DNA_repair")

dif_expr_BPSC=function(m, l.bp4, cond){
  
  ind=grep(cond,colnames(m))
  m1=m[,ind] ; m2=m[,-c(ind)]
  
  ## Remplir les NA
  tmp1 <- m1
  for (i in 1:nrow(tmp1)) {
    n <- length(which(is.na(tmp1[i,])))
    if(is.na(l.bp4[[1]][i,1])==FALSE){
      if (length(which(is.na(tmp1[i,]))) > 0){
        tmp1[i,which(is.na(tmp1[i,]))] <- rBP(n, l.bp4[[1]][i,])
      }
    }else{
      tmp1[i,which(is.na(tmp1[i,]))]<-NA
    }
  }
  
  tmp2 <- m2
  for (i in 1:nrow(tmp2)) {
    n <- length(which(is.na(tmp2[i,])))
    if(is.na(l.bp4[[2]][i,1])==FALSE){
      if (length(which(is.na(tmp2[i,]))) > 0){
        tmp2[i,which(is.na(tmp2[i,]))] <- rBP(n, l.bp4[[2]][i,])
      }
    }else{
      tmp2[i,which(is.na(tmp2[i,]))]<-NA
    }
  }
  
  tmp <- cbind(tmp1, tmp2)
  ctidx <- rep(c(1,0), c(ncol(m1),ncol(m2)))
  design=model.matrix(~ctidx)
  res=BPglm(data=tmp, design=design, coef=2, estIntPar=TRUE)
  
  matrice_diff=matrix(ncol=1,nrow=nrow(m))
  rownames(matrice_diff)=rownames(m)
  colnames(matrice_diff)=c(paste(cond,"_Control",sep=""))
  
  for(i in 1:length(res$PVAL)){
    if((res$PVAL[i]<0.01)==T & is.na(res$PVAL[i])==F){
      
      matrice_diff[i,1]=ifelse(mean(m1[i,],na.rm=T)>mean(m2[i,],na.rm=T),1,-1)
      
    }
    else(matrice_diff[i,1]=ifelse(is.na(res$PVAL[i])==T,"erreur",0))
  }
  
  p.value=res$PVAL
  l=list("p-value"=p.value,"Matrice"=matrice_diff)
  return(l)
}

diff_glyco=dif_expr_BPSC(data.glyco.thresh,coeff_BP_glyco.thresh,"Glyco")
diff_EMT=dif_expr_BPSC(data.EMT.thresh,coeff_BP_EMT.thresh,"EMT")
diff_DNA_repair=dif_expr_BPSC(data.DNA_repair.thresh,coeff_BP_DNA_repair.thresh,"DNA_repair")


gene.dif.glyco.thresh=rownames(diff_glyco[[2]])[which(diff_glyco[[2]]!=0)]
gene.dif.EMT.thresh=rownames(diff_EMT[[2]])[which(diff_EMT[[2]]!=0)]
gene.dif.DNA_repair.thresh=rownames(diff_DNA_repair[[2]])[which(diff_DNA_repair[[2]]!=0)]


###################################       Combinaisons    #######################################


combinaison=function(liste,nbmax){
  l=list()
  for(i in 2:min(nbmax,length(liste))){
    l[[c(paste("Combinaison de",i,"gènes"))]]=combn(liste,i)
  }
  return(l)
}

comb_glyco.2=combinaison(gene.dif.glyco.thresh,10)
comb_EMT.2=combinaison(gene.dif.EMT.thresh,10)
comb_DNA_repair.2=combinaison(gene.dif.DNA_repair.thresh,10)


m.nom=function(comb){
  v=list()
  for(j in 1:length(comb)){
    v[[j]]=matrix(ncol=1,nrow=ncol(comb[[j]]))
    for(i in 1:ncol(comb[[j]])){
      x=comb[[j]][1,i]
      for(k in 2:nrow(comb[[j]])){
        x=str_c(x,comb[[j]][k,i],sep=", ")
      }
      v[[j]][i,]=x
    }
  }
  names(v)=names(comb)
  return(v)
}

n.comb_glyco.2=m.nom(comb_glyco.2)
n.comb_EMT.2=m.nom(comb_EMT.2)
n.comb_DNA_repair.2=m.nom(comb_DNA_repair.2)



list.mcomb=function(proba,comb,ncomb,cond){
  l=list() ; l1=list() ; l2=list()
  
  for(i in 1:length(comb)){
    k=nrow(comb[i])
    l1[[i]]=list() ; l2[[i]]=list() 
    
    for(j in 1:ncol(comb[[i]])){
      n.gene=comb[[i]][,j]
      l1[[i]][[j]]=as.matrix(c(proba[[n.gene[1]]][,6]))
      l2[[i]][[j]]=as.matrix(c(proba[[n.gene[1]]][,7]))
      
      for(z in 2:length(n.gene)){
        l1[[i]][[j]]=cbind(l1[[i]][[j]],c(proba[[n.gene[z]]][,6]))
        l2[[i]][[j]]=cbind(l2[[i]][[j]],c(proba[[n.gene[z]]][,7]))
      }
      colnames(l1[[i]][[j]])=c(n.gene) ; colnames(l2[[i]][[j]])=c(n.gene)
      
    }
    names(l1[[i]])=ncomb[[i]] ; names(l2[[i]])=ncomb[[i]]
  }
  l=list(l1,l2)
  names(l)=c(cond,"Control")
  names(l[[1]])=names(comb) ; names(l[[2]])=names(comb)
  return(l)
}

l.mcomb_glyco.2=list.mcomb(likelihood_glyco.thresh,comb_glyco.2,n.comb_glyco.2,"Glycolysis")
l.mcomb_EMT.2=list.mcomb(likelihood_EMT.thresh,comb_EMT.2,n.comb_EMT.2,"EMT")
l.mcomb_DNA_repair.2=list.mcomb(likelihood_DNA_repair.thresh,comb_DNA_repair.2,
                                n.comb_DNA_repair.2,"DNA_repair")
coeff_glm=function(lmcomb,cond){
  m1=list(); m2=list()
  glm1=list() ; glm2=list()
  
  l1=lmcomb[[1]] ; l2=lmcomb[[2]]
  
  n1=length(grep(cond,rownames(l1[[1]][[1]])))
  n2=nrow(l1[[1]][[1]])-n1
  
  
  for(i in 1:length(l1)){
    print(names(l1)[i])
    m1[[i]]=list(); m2[[i]]=list()
    glm1[[i]]=list() ; glm2[[i]]=list()
    for(j in 1:length(l1[[i]])){
      m1[[i]][[j]]=list() ;  m2[[i]][[j]]=list()
      glm1[[i]][[j]]=list() ;  glm2[[i]][[j]]=list()
      
      x1=as.data.frame(l1[[i]][[j]])
      x2=as.data.frame(l2[[i]][[j]])
      resp1=c(rep(1,n1),rep(0,n2))
      resp2=c(rep(0,n1),rep(1,n2))
      
      m1[[i]][[j]]=matrix(nrow=nrow(x1),ncol=(ncol(x1)+1))
      m2[[i]][[j]]=matrix(nrow=nrow(x2),ncol=(ncol(x2)+1))
      
      colnames(m1[[i]][[j]])=c("Intercept",colnames(l1[[i]][[j]]))
      rownames(m1[[i]][[j]])=c(rownames(l1[[i]][[j]]))
      
      colnames(m2[[i]][[j]])=c("Intercept",colnames(l2[[i]][[j]]))
      rownames(m2[[i]][[j]])=c(rownames(l2[[i]][[j]]))
      
      for(k in 1:nrow(x1)){
        glm1[[i]][[j]][[k]]=list() ;  glm2[[i]][[j]][[k]]=list()
        
        y1=x1[-k, ]
        glm1[[i]][[j]][[k]]=glm(resp1[-k]~.,family="binomial",data=y1)
        m1[[i]][[j]][k,]=c(glm1[[i]][[j]][[k]]$coefficients)
        
        
        y2=x2[-k, ]
        glm2[[i]][[j]][[k]]=glm(resp2[-k]~.,family="binomial",data=y2)
        m2[[i]][[j]][k,]=c(glm2[[i]][[j]][[k]]$coefficients)
      }
      names(glm1[[i]][[j]])=names(glm2[[i]][[j]])=rownames(l1[[i]][[j]])
    }
    names(m1[[i]])=names(l1[[i]])
    names(m2[[i]])=names(l2[[i]])
    names(glm1[[i]])=names(l1[[i]])
    names(glm2[[i]])=names(l2[[i]])
  }
  names(m1)=names(l1)
  names(m2)=names(l2)
  names(glm1)=names(l1)
  names(glm2)=names(l2)
  
  coeff=list(m1,m2)
  l.glm=list(glm1,glm2)
  names(coeff)=names(lmcomb)
  names(l.glm)=names(lmcomb)
  
  v=list(coeff,l.glm)
  names(v)=c("Coefficient","Modèle")
  return(v)
}

coeff_glm_glyco.thresh.2=coeff_glm(l.mcomb_glyco.2,"Glyco")
coeff_glm_EMT.thresh.2=coeff_glm(l.mcomb_EMT.2,"EMT")
coeff_glm_DNA_repair.thresh.2=coeff_glm(l.mcomb_DNA_repair.2,"DNA_repair")


p.logit=function(proba,model){
  l1=proba[[1]] ; l2=proba[[2]]
  
  logit1=list() ; logit2=list()
  
  for(i in 1:length(l1)){
    logit1[[i]]=list() ; logit2[[i]]=list()
    for(j in 1:length(l1[[i]])){
      logit1[[i]][[j]]=vector("numeric",nrow(l1[[i]][[j]]))
      logit2[[i]][[j]]=vector("numeric",nrow(l2[[i]][[j]]))
      
      for(k in 1:nrow(l1[[i]][[j]])){
        x=as.data.frame(t(l1[[i]][[j]][k,]))
        y=as.data.frame(t(l2[[i]][[j]][k,]))
        
        logit1[[i]][[j]][k]=predict(model[[1]][[i]][[j]][[k]],x)
        logit2[[i]][[j]][k]=predict(model[[2]][[i]][[j]][[k]],y)
      }
      names(logit1[[i]][[j]])=names(logit2[[i]][[j]])=rownames(l1[[i]][[j]])
    }
    names(logit1[[i]])=names(logit2[[i]])=names(l1[[i]])
    
  }
  names(logit1)=names(logit2)=names(l1)
  
  l=list(logit1,logit2)
  names(l)=names(proba)
  return(l)
}

glm.glyco.thresh.2=p.logit(l.mcomb_glyco.2,coeff_glm_glyco.thresh.2[[2]])
glm.EMT.thresh.2=p.logit(l.mcomb_EMT.2,coeff_glm_EMT.thresh.2[[2]])
glm.DNA_repair.thresh.2=p.logit(l.mcomb_DNA_repair.2,coeff_glm_DNA_repair.thresh.2[[2]])


antilogit=function(x){
  return(exp(x)/(1+exp(x)))
}

p.glm=function(glm){
  l1=glm[[1]] ; l2=glm[[2]]
  m1=list() ; m2=list()
  for(i in 1:length(l1)){
    
    m=matrix(ncol=length(l1[[i]]),nrow=length(l1[[i]][[1]]))
    
    for(j in 1:length(l1[[i]])){
      m[,j]=sapply(l1[[i]][[j]],antilogit)
    }
    m1[[i]]=m
    colnames(m1[[i]])=names(l1[[i]])
    rownames(m1[[i]])=names(l1[[i]][[1]])
  }
  
  for(i in 1:length(l2)){
    m=matrix(ncol=length(l2[[i]]),nrow=length(l2[[i]][[1]]))
    for(j in 1:length(l2[[i]])){
      m[,j]=sapply(l2[[i]][[j]],antilogit)
    }
    m2[[i]]=m
    colnames(m2[[i]])=names(l2[[i]])
    rownames(m2[[i]])=names(l2[[i]][[1]])
    
  }
  names(m1)=names(m2)=names(l1)
  l=list(m1,m2)
  names(l)=names(glm)
  return(l)
}

proba.glyco.glm.thresh.2=p.glm(glm.glyco.thresh.2)
proba.EMT.glm.thresh.2=p.glm(glm.EMT.thresh.2)
proba.DNA_repair.glm.thresh.2=p.glm(glm.DNA_repair.thresh.2)


### Courbes ROC
roc.glm=function(proba,cond){
  l1=proba[[1]] ; l2=proba[[2]]
  cond1=names(proba)[1]
  cond2=names(proba)[2]
  
  roc1=list() ; roc2=list()
  
  for(i in 1:length(l1)){
    pdf(paste("ROC_glm_",names(l1)[i],cond1,".pdf",sep=""))
    
    n1=length(grep(cond,rownames(l1[[i]])))
    n2=nrow(l1[[i]])-n1
    
    roc1[[i]]=list() ; roc2[[i]]=list()
    
    for(j in 1:ncol(l1[[i]])){
      par(mar=c(7,5,3,3),xpd=F)
      roc1[[i]][[j]]=list() ; roc2[[i]][[j]]=list()
      
      roc1[[i]][[j]]=roc(response= c(rep(1,n1),rep(0,n2)), predictor=l1[[i]][,j])
      roc2[[i]][[j]]=roc(response= c(rep(0,n1),rep(1,n2)), predictor=l2[[i]][,j])
      plot(1-roc1[[i]][[j]]$specificities, roc1[[i]][[j]]$sensitivities, type="l", 
           xlab="1 - Specificity", ylab="Sensitivity", col="blue", 
           main=paste(cond1,"from",cond2,sep=" "),ylim=c(0,1.05))
      points(1-roc2[[i]][[j]]$specificities, roc2[[i]][[j]]$sensitivities, type="l", col="red")
      points(x=c(0,1), y=c(0,1), type="l", lty=2)
      legend("bottomright",inset=0.07, col=c("blue", "red"), legend=c(cond1, cond2), lty=1,
             box.lty = 0,cex=0.9,horiz=T)
      
      legend("top",inset=0.01,legend=colnames(l1[[i]])[j],
             text.font=4,box.lty=0,cex=0.6)
      
      par(xpd=T)
      text(0.5,-0.25,paste("AUC",cond1,"=", round(auc(roc1[[i]][[j]]),digits = 3),"   ", 
                           "AUC",cond2, "=",round(auc(roc2[[i]][[j]]),digits = 3)),font=2, cex=0.8)
      
      
    }
    names(roc1[[i]])=names(roc2[[i]])=colnames(l1[[i]])
    dev.off()
  }
  names(roc1)=names(roc2)=names(l1)
  l=list(roc1,roc2)
  names(l)=c(cond1,cond2)
  return(l)
}

roc_glm_glyco.thresh.2=roc.glm(proba.glyco.glm.thresh.2,"Glyco")
roc_glm_EMT.thresh.2=roc.glm(proba.EMT.glm.thresh.2,"EMT")
roc_glm_DNA_repair.thresh.2=roc.glm(proba.DNA_repair.glm.thresh.2,"DNA_repair")



### Récupération des aires sous la courbe : AUC
calcul.auc.1=function(roc){
  l=list()
  roc1=roc[[1]] ; roc2=roc[[2]]
  for(i in 1:length(roc1)){
    m.auc=matrix(ncol=2,nrow=length(roc1[[i]]))
    rownames(m.auc)=names(roc1[[i]])
    colnames(m.auc)=names(roc)
    
    for(j in 1:length(roc1[[i]])){
      m.auc[j,]=c(roc1[[i]][[j]]$auc,roc2[[i]][[j]]$auc)
    }
    l[[i]]=m.auc
  }
  names(l)=names(roc1)
  
  mat=NULL
  for(i in 1:length(l)){
    mat=rbind(mat,l[[i]])
  }
  
  return(mat)
}

auc_glm_glyco.thresh.2=calcul.auc.1(roc_glm_glyco.thresh.2)
auc_glm_EMT.thresh.2=calcul.auc.1(roc_glm_EMT.thresh.2)
auc_glm_DNA_repair.thresh.2=calcul.auc.1(roc_glm_DNA_repair.thresh.2)


ALL_auc_glm_glyco.thresh.2=rbind(auc_glyco.thresh[which(diff_glyco[[2]]!=0),],
                                 auc_glm_glyco.thresh.2)
ALL_auc_glm_EMT.thresh.2=rbind(auc_EMT.thresh[which(diff_EMT[[2]]!=0),],
                               auc_glm_EMT.thresh.2)
ALL_auc_glm_DNA_repair.thresh.2=rbind(auc_DNA_repair.thresh[which(diff_DNA_repair[[2]]!=0),],
                                      auc_glm_DNA_repair.thresh.2)


vrep.glyco.thresh.2=c(rep(1,4),rep(2,6),rep(3,4),rep(4,1))
vrep.EMT.thresh.2=c(rep(1,3),rep(2,3),rep(3,1))
vrep.DNA_repair.thresh.2=c(rep(1,3),rep(2,3),rep(3,1))

## Evolution AUC en fonction du nombre de gènes
library(beeswarm)
evolution.auc=function(auc,vrep,chaine,cond,ind,ind2){
  pdf(paste("Evolution_AUC_",cond,"_",chaine,".pdf",sep=""))
  par(mar=c(4.1,5.1,5,4.1),xpd=T)
  x=cbind(auc[,1],vrep)
  colnames(x)=c("X1","X2")
  boxplot(X1~X2, data=x,xlab="Combinaison de i gènes",ylab="AUC",
          ylim=c(min(x[,1],na.rm=T),max(x[,1],na.rm=T)),
          main="Evolution des AUC en fonction du nombre de gènes",outline=F,
          col=rainbow(10,alpha=0.3))
  beeswarm(X1~X2,add=T, col=rainbow(9), pch=20, corral="wrap",data=x)
  
  text(ind,max(auc[,1],na.rm=T)+ind2,chaine,cex=0.8,font=2)
  dev.off()
}

evolution.auc(ALL_auc_glm_glyco.thresh.2,vrep.glyco.thresh.2,"Régression logistique","Glycolysis",2.5,0.02)
evolution.auc(ALL_auc_glm_EMT.thresh.2,vrep.EMT.thresh.2,"Régression logistique","EMT",2,0.01)
evolution.auc(ALL_auc_glm_DNA_repair.thresh.2,vrep.DNA_repair.thresh.2,"Régression logistique","DNA_repair",2,0.005)



### Prédictions

m_AUC_glm_glyco.thresh=cbind(ALL_auc_glm_glyco.thresh.2[,1],vrep.glyco.thresh.2)
m_AUC_glm_EMT.thresh=cbind(ALL_auc_glm_EMT.thresh.2[,1],vrep.EMT.thresh.2)
m_AUC_glm_DNA_repair.thresh=cbind(ALL_auc_glm_DNA_repair.thresh.2[,1],vrep.DNA_repair.thresh.2)



find.AUCmax=function(m,x){
  ind=which(m[,2]==x)
  v=m[ind,1]
  return(v[which(v==max(v))])
}

#combMAX1_glyco.thresh=find.AUCmax(m_AUC_glm_glyco.thresh,1)
combMAX2_glyco.thresh=find.AUCmax(m_AUC_glm_glyco.thresh,2)
combMAX3_glyco.thresh=find.AUCmax(m_AUC_glm_glyco.thresh,3)
combMAX4_glyco.thresh=m_AUC_glm_glyco.thresh[15,1]
names(combMAX4_glyco.thresh)=rownames(m_AUC_glm_glyco.thresh)[15]

combMAX1_EMT.thresh=find.AUCmax(m_AUC_glm_EMT.thresh,2)
combMAX2_EMT.thresh=find.AUCmax(m_AUC_glm_EMT.thresh,2)
combMAX3_EMT.thresh=find.AUCmax(m_AUC_glm_EMT.thresh,3)
names(combMAX3_EMT.thresh)=rownames(m_AUC_glm_EMT.thresh)[7]

combMAX2_DNA_repair.thresh=find.AUCmax(m_AUC_glm_DNA_repair.thresh,2)
combMAX3_DNA_repair.thresh=find.AUCmax(m_AUC_glm_DNA_repair.thresh,3)
names(combMAX3_DNA_repair.thresh)=rownames(m_AUC_glm_DNA_repair.thresh)[7]


return_proba=function(proba,comb,x){
  ind=grep(names(comb),colnames(proba[[1]][[x-1]]))
  return(proba[[1]][[x-1]][,ind])
}

## Avec plusieurs AUC max
return_proba_2=function(proba,comb,x){
  v=list()
  ind=NULL
  for(i in 1:length(comb)){
    ind[i]=grep(names(comb)[i],colnames(proba[[1]][[x-1]]))
    v[[i]]=proba[[1]][[x-1]][,ind[i]]
  }
  names(v)=names(comb)
  return(v)
}


p_glm2_glyco.thresh=return_proba(proba.glyco.glm.thresh.2,combMAX2_glyco.thresh,2)
p_glm3_glyco.thresh=return_proba(proba.glyco.glm.thresh.2,combMAX3_glyco.thresh,3)
p_glm4_glyco.thresh=return_proba(proba.glyco.glm.thresh.2,combMAX4_glyco.thresh,4)

p_glm2_EMT.thresh=return_proba_2(proba.EMT.glm.thresh.2,combMAX2_EMT.thresh,2)
p_glm3_EMT.thresh=return_proba(proba.EMT.glm.thresh.2,combMAX3_EMT.thresh,3)

p_glm2_DNA_repair.thresh=return_proba(proba.DNA_repair.glm.thresh.2,combMAX2_DNA_repair.thresh,2)
p_glm3_DNA_repair.thresh=return_proba(proba.DNA_repair.glm.thresh.2,combMAX3_DNA_repair.thresh,3)


l_pred_glyco.thresh=list(p_glm2_glyco.thresh,p_glm3_glyco.thresh,p_glm4_glyco.thresh)
names(l_pred_glyco.thresh)=c(names(combMAX2_glyco.thresh),names(combMAX3_glyco.thresh),
                             names(combMAX4_glyco.thresh))

l_pred_EMT.thresh=list(p_glm2_EMT.thresh[[1]],p_glm2_EMT.thresh[[2]],p_glm3_EMT.thresh)
names(l_pred_EMT.thresh)=c(names(combMAX2_EMT.thresh),names(combMAX3_EMT.thresh))

l_pred_DNA_repair.thresh=list(p_glm2_DNA_repair.thresh,p_glm3_DNA_repair.thresh)
names(l_pred_DNA_repair.thresh)=c(names(combMAX2_DNA_repair.thresh),names(combMAX3_DNA_repair.thresh))



plot_prediction=function(l,cond, title=NULL, log=T){
  if (is.null(title))
    title=paste("prediction_glm_", cond,"_log10.pdf",sep="")
  pdf(title)
  par(cex.axis=0.6,font.axis=2,mar=c(5,5,5,5),mfrow=c(3,1))
  
  for(i in 1:length(l)){
    m <- l[[i]]
    if (log==T)
      m=log10(m)
    image(as.matrix(m),col=rev(grey.colors(100)),axes=F, frame.plot=TRUE,ylab=names(l)[i],font.lab=2)
    axis(1, seq(0,1,1/35),names(l[[i]]),las=2)
  }
  dev.off()
}

plot_prediction(l_pred_glyco.thresh,"Glycolysis")
plot_prediction(l_pred_EMT.thresh,"EMT")
plot_prediction(l_pred_DNA_repair.thresh,"DNA_repair")

l_pred_EMT.all=list(p_glm2_EMT.thresh,p_glm2_EMT.thresh,p_glm3_EMT.thresh)
plot_prediction(l_pred_EMT.thresh,"EMT", "prediction_glm_EMT_all_log10.pdf", log=T)

### Sensibilité et spécificité
calcul_sensi_speci=function(l,seuil,cond){
  ind=grep(cond,names(l[[1]]))
  n1=length(ind)
  n2=length(l[[1]])-n1
  
  tab=matrix(ncol=2,nrow=length(l))
  colnames(tab)=c("Sensibilité","Spécificité")
  rownames(tab)=names(l)
  
  mat=list()
  
  for(i in 1:length(l)){
    mat[[i]]=matrix(ncol=3,nrow=3)
    rownames(mat[[i]])=c("+","-","TOTAL")
    colnames(mat[[i]])=c(cond,"Control","TOTAL")
    
    x=which(l[[i]]>seuil)
    
    mat[[i]][1,3]=length(x)
    mat[[i]][3,3]=length(l[[i]])
    mat[[i]][3,c(1,2)]=c(n1,n2)
    mat[[i]][2,3]=mat[[i]][3,3]-mat[[i]][1,3]
    
    x=which(l[[i]]>seuil)
    
    n3=length(grep(cond,names(x)))
    mat[[i]][1,1]=n3
    
    mat[[i]][2,1]=mat[[i]][3,1]-mat[[i]][1,1]
    mat[[i]][1,2]=mat[[i]][1,3]-mat[[i]][1,1]
    mat[[i]][2,2]=mat[[i]][2,3]-mat[[i]][2,1]
    
    tab[i,]=c((mat[[i]][1,1]/mat[[i]][3,1]),(mat[[i]][2,2]/mat[[i]][3,2]))
  }
  
  names(mat)=names(l)
  
  y=list(mat,tab)
  names(y)=c("Tableau","Sensibilité-Spécificité")
  return(y)
  
}

glyco.ss = calcul_sensi_speci(l_pred_glyco.thresh,0.1,"Glyco")[["Sensibilité-Spécificité"]]
glyco.final = cbind(m_AUC_glm_glyco.thresh[rownames(glyco.ss),1,drop=F], glyco.ss)
colnames(glyco.final)[1]="AUC"
write.table(glyco.final, file="glycolysis_AUC_sens_spec.txt", sep="\t", quote=F, col.names=NA)

EMT.ss = calcul_sensi_speci(l_pred_EMT.thresh,0.1,"EMT")[["Sensibilité-Spécificité"]]
EMT.final = cbind(m_AUC_glm_EMT.thresh[rownames(EMT.ss),1,drop=F], EMT.ss)
colnames(EMT.final)[1]="AUC"
write.table(EMT.final, file="EMT_AUC_sens_spec.txt", sep="\t", quote=F, col.names=NA)

DNArepair.ss = calcul_sensi_speci(l_pred_DNA_repair.thresh,0.1,"DNA_repair")[["Sensibilité-Spécificité"]]
DNArepair.final = cbind(m_AUC_glm_DNA_repair.thresh[rownames(DNArepair.ss),1,drop=F], DNArepair.ss)
colnames(DNArepair.final)[1]="AUC"
write.table(DNArepair.final, file="DNArepair_AUC_sens_spec.txt", sep="\t", quote=F, col.names=NA)
