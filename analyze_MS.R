# CE = Cholesteryl Esters
# LPC= Lysophosphocholine
# PA= phosphatic acid
# PS= phosphatidylglycerol
# PC= phosphocholine
# PE= Phospoethanolamine
# SM= sphyngomylin
# Cer= Ceramides
# DG= Diacylglycerol
# TG=Triacylglycerol

msFreeArachidonicAcid <- "[FFA(20:4)+OAc]"

library(Rtsne)
library(limma)
library(xlsx)
#install.packages("Rtsne")


######################################################################
############# Read the MS data #######################################
######################################################################

## Read condition matrix
mscond <- read.csv("ms/mscond.csv",stringsAsFactors = FALSE)
rownames(mscond) <- mscond$sample
colisko <- rep("black",nrow(mscond))
colisko[mscond$isko==1] <- "red"

## Read positive
msdatapos <- read.xlsx("ms/T cells Pos Lipid.xlsx",sheetIndex = 2,stringsAsFactors = FALSE,header = FALSE)#[-(1:2),]
msdataposnames <- msdatapos[1,-1]
msdataposmz <- msdatapos[2,-1]
rownames(msdatapos)<-msdatapos$X1
msdatapos <- msdatapos[-(1:3),-1]
colnames(msdatapos) <- msdataposnames

## Read negative
msdataneg <- read.xlsx("ms/T cells Neg Lipid.xlsx",sheetIndex = 2,stringsAsFactors = FALSE,header = FALSE)#[-(1:2),]
msdatanegnames <- msdataneg[1,-1]
msdatanegmz <- msdataneg[2,-1]
rownames(msdataneg)<-msdataneg$X1
msdataneg <- msdataneg[-(1:3),-1]
colnames(msdataneg) <- msdatanegnames

## Combine positive and negative
msdataispos <- c(rep(FALSE,length(msdatanegnames)), rep(TRUE,length(msdatanegnames)))
msdatanames <- cbind(msdatanegnames, msdataposnames)
msdatamz <- unlist(c(msdatanegmz, msdataposmz))
class(msdatamz) <- "double"
msdata <- cbind(msdataneg, msdatapos)
colnames(msdata) <- 1:ncol(msdata)
msdata <- as.matrix(msdata)
class(msdata)<-"double"


#size-factor normalization
sfn <- estimateSizeFactorsForMatrix(t(msdata))
for(i in 1:nrow(msdata)){
  msdata[i,] <- msdata[i,]/sfn[i]
}

######################################################################
############## Perform PCA ###########################################
######################################################################


doMSpca <- function(p1=1, p2=2){
  mod<-princomp(t(msdata[keep,]), cor=FALSE)
  par(mfrow=c(1,1))
  plot(mod$loadings[,p1],mod$loadings[,p2],cex=0)
  text(label=mscond$sample[keep],mod$loadings[,p1],mod$loadings[,p2], col=colisko[keep])
}

pdf("ms/out_PCA_all.pdf")
keep <- mscond$mousenum>0 #Clusters by: in vivo / late time / early time. then by individual
doMSpca(3,4)
dev.off()

pdf("ms/out_PCA_invivo.pdf")
keep <- mscond$mousenum>5 #Clusters by individual
doMSpca(3,4)
dev.off()

pdf("ms/out_PCA_invitro.pdf")
keep <- mscond$mousenum<5 #Clusters by: time, then possibly individual
doMSpca(3,4)
dev.off()

######################################################################
############## Perform tSNE ##########################################
######################################################################


pdf("ms/out_tsne.pdf")
set.seed(123) # set random seed
rtsne_out <- Rtsne(as.matrix(msdata), pca = FALSE, verbose = TRUE, perplexity = 7)
plot(rtsne_out$Y, asp = 1, pch = 20, col = "blue", 
     cex = 0, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, 
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", 
     main = "2D t-SNE projection")
text(rtsne_out$Y, labels = mscond$sample, col=colisko)
dev.off()

# Here it clusters on individual, then *maybe* in KO-status. but equally well on tissue.
# Essentially no difference between in vitro Th types. But plenty difference with time




######################################################################
######## Compare groups by straight linear models ####################
######################################################################


simptest <- function(keep, form, number=10, coef=NULL){
  mscond <- mscond[keep,]
  mscond$mousenum <- factor(mscond$mousenum)
  mscond$time <- factor(mscond$time)
  mscond$tissue <- factor(mscond$tissue)
  
  phenoData <- new("AnnotatedDataFrame", data=mscond)
  #msdata <- as.matrix(msdata)
  assdata <- t(as.matrix(msdata[keep,]))
  eset <- ExpressionSet(assdata, phenoData = phenoData)
  design <- model.matrix(form, pData(eset))
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)  
  vt <- topTable(fit, number=ncol(msdata), sort.by = "none", coef = coef)  
  vt$m.z <- msdatamz
  
  vt <- as.data.frame(vt)
  names(msdatanames) <- NULL
  vt <- data.frame(vt, name=unlist(msdatanames[as.integer(rownames(vt))]))
  vt  
}

writeMummichog <- function(st, ispos, fname="ms/out_list.csv"){
  fname <- "ms/out_list.csv"
  ispos <- TRUE
  i <- st$P.Value<0.05 & ispos==msdataispos
  sum(i)
  
  d <- data.frame(m.z=unlist(msdatamz[i]), p.value=st$P.Value[i])
  write.table(file=fname, d, row.names = FALSE, quote = FALSE, sep="\t")
}

#look for FFA 20:4

allLinMod <- function(){
  mscond$mousenum
  msdata
  
}




st <- simptest(mscond$tissue=="spleen", ~1+isko, coef = "isko")[,c("name","m.z","logFC","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:10],]
st[st$name==msFreeArachidonicAcid,]

st <- simptest(mscond$tissue=="lung", ~1+isko, coef="isko")[,c("name","m.z","logFC","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:10],]


st <- simptest(mscond$tissue %in% c("spleen","lung"), ~1+tissue+isko, coef="isko")[,c("name","m.z","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:20],]


st <- simptest(mscond$istc==1, ~1+mousenum+time,number=30,coef="time")[,c("name","m.z","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:20],]



st <- simptest(mscond$istc==1, ~1+time, coef="time")[,c("name","m.z","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:20],]


st <- simptest(mscond$time==148 & mscond$tissue %in% c("th0","th2"), ~1+mousenum+tissue, coef="tissue")[,c("name","m.z","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:20],]



st <- simptest(mscond$time==148 & mscond$tissue %in% c("th1","th2"), ~1+mousenum+tissue, coef="tissue")[,c("name","m.z","P.Value","adj.P.Val")]
st[order(st$P.Value)[1:20],]


#st <- simptest(mscond$istc==1, ~1+time, coef="tissue")[,c("name","m.z","P.Value","adj.P.Val")]
st_th02 <- simptest(mscond$time==148 & mscond$tissue %in% c("th0","th2"), ~tissue+mousenum)[,c("name","m.z","P.Value","adj.P.Val")]
st_th01 <- simptest(mscond$time==148 & mscond$tissue %in% c("th0","th1"), ~tissue+mousenum)[,c("name","m.z","P.Value","adj.P.Val")]
st_th12 <- simptest(mscond$time==148 & mscond$tissue %in% c("th1","th2"), ~tissue+mousenum)[,c("name","m.z","P.Value","adj.P.Val")]


mscond$tissue

######################################## 
######################################## 
######################################## later, in addition to eBayes
######################################## 
######################################## 
######################################## 
######################################## 
fit <- lmFit(M)

#  Ordinary t-statistic
par(mfrow=c(1,2))
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma




#########################################################################
###################### Compare KO vs WT, in vivo ########################
#########################################################################


keepMetabolite <- which(apply(msdata,2,mean)>=0)
length(keepMetabolite)

keep <- mscond$tissue %in% c("spleen","lung")
form <- ~1+tissue+isko

### Retain some samples
mscond_red <- mscond[keep,]
mscond_red$mousenum <- factor(mscond_red$mousenum)
mscond_red$time <- factor(mscond_red$time)
mscond_red$tissue <- factor(mscond_red$tissue)


phenoData <- new("AnnotatedDataFrame", data=mscond_red)
assdata <- t(as.matrix(msdata[keep,keepMetabolite]))
eset <- ExpressionSet(assdata, phenoData = phenoData)
design <- model.matrix(form, pData(eset))


fit <- lmFit(eset, design)
fit <- eBayes(fit)  
vt <- topTable(fit, number=ncol(msdata[,keepMetabolite]), sort.by = "none", coef = "isko")  
vt$m.z <- msdatamz[keepMetabolite]

vt <- as.data.frame(vt)
names(msdatanames) <- NULL
vt$name <- (unlist(msdatanames))[keepMetabolite]
#vt <- data.frame(vt, name=msdatanames[keepMetabolite])
#vt <- data.frame(vt, name=unlist(msdatanames[as.integer(rownames(vt))]))
vt <- vt[order(vt$P.Value),]
vt[abs(vt$logFC)>0.01,]  

msFreeArachidonicAcid
vt[vt$name==msFreeArachidonicAcid,]
vt[grep("FFA",vt$name),]
vt

  







#########################################################################
###################### Compare Th1 vs Th2, in vitro ########################
#########################################################################

# 
# keepMetabolite <- which(apply(msdata,2,mean)>=0)
# length(keepMetabolite)
# 
# keep <- mscond$tissue %in% c("th2","th1")
# form <- ~1+mousenum+tissue
# 
# ### Retain some samples
# mscond_red <- mscond[keep,]
# mscond_red$mousenum <- factor(mscond_red$mousenum)
# mscond_red$time <- factor(mscond_red$time)
# mscond_red$tissue <- factor(mscond_red$tissue)
# 
# 
# phenoData <- new("AnnotatedDataFrame", data=mscond_red)
# assdata <- t(as.matrix(msdata[keep,keepMetabolite]))
# eset <- ExpressionSet(assdata, phenoData = phenoData)
# design <- model.matrix(form, pData(eset))
# 
# 
# fit <- lmFit(eset, design)
# fit <- eBayes(fit)  
# vt <- topTable(fit, number=ncol(msdata[,keepMetabolite]), sort.by = "none", coef = "tissueth2")  
# vt$m.z <- msdatamz[keepMetabolite]
# 
# vt <- as.data.frame(vt)
# names(msdatanames) <- NULL
# vt$name <- (unlist(msdatanames))[keepMetabolite]
# vt <- vt[order(vt$P.Value),]
# vt



# Th1 vs Th2, with mouse batch
# logFC      AveExpr            t     P.Value adj.P.Val         B       m.z                                               name
# 17  -6.579263e-02 2.611265e-01 -7.033071344 0.008812544 0.7021615 -2.019808  904.7030                       [PC(40:0)+OAc][PE(43:0)+OAc]
# 205 -4.039712e-02 5.031960e-02 -6.695510373 0.009982287 0.7021615 -2.176500  922.7853                                     [TG(56:7)+NH4]
# 204 -5.130124e-02 6.931651e-02 -5.886515150 0.013797252 0.7021615 -2.587284  924.8005                                     [TG(56:6)+NH4]
# 201 -5.902283e-02 7.821056e-02 -5.816118752 0.014217723 0.7021615 -2.625624  926.8153                                     [TG(56:5)+NH4]
# 168 -8.946336e-02 4.263329e-01 -5.448851580 0.016718477 0.7021615 -2.833114  848.7712                                     [TG(50:2)+NH4]
# 166 -1.238143e-01 4.452295e-01 -4.893212640 0.021773107 0.7021615 -3.173003  850.7850                                     [TG(50:1)+NH4]
# 214 -3.186529e-02 3.043348e-02 -4.823551189 0.022546470 0.7021615 -3.218021  948.8027                                     [TG(58:8)+NH4]




#########################################################################
###################### Compare Th1 vs Th2, in vitro ########################
#########################################################################


keepMetabolite <- which(apply(msdata,2,mean)>=0)
length(keepMetabolite)

doMStest <- function(keep, form, mscond_orig=mscond){
  ### Retain some samples
  mscond_red <- mscond_orig[keep,]
  mscond_red$mousenum <- factor(mscond_red$mousenum)
  mscond_red$time <- factor(mscond_red$time)
  mscond_red$tissue <- factor(mscond_red$tissue)
  
  
  phenoData <- new("AnnotatedDataFrame", data=mscond_red)
  assdata <- t(as.matrix(msdata[keep,keepMetabolite]))
  eset <- ExpressionSet(assdata, phenoData = phenoData)
  design <- model.matrix(form, pData(eset))
  
  
  
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)  
  vt <- topTable(fit, number=ncol(msdata[,keepMetabolite]), sort.by = "none", coef = rev(colnames(fit$coefficients))[1])  
  vt$m.z <- msdatamz[keepMetabolite]
  
  vt <- as.data.frame(vt)
  names(msdatanames) <- NULL
  vt$name <- (unlist(msdatanames))[keepMetabolite]
  vt <- vt[order(vt$P.Value),]
  vt
}


keep <- mscond$tissue %in% c("th2","th0")
form <- ~1+mousenum+tissue
doMStest(keep, form)

# logFC      AveExpr            t    P.Value adj.P.Val         B       m.z                                               name
# 73   8.216168e-02 8.659120e-02  2.804354658 0.08380719 0.6601969 -4.687132  720.5555                           [PC(31:0)+H][PE(34:0)+H]
# 8    7.904578e-01 1.174005e+00  2.702556387 0.09024032 0.6601969 -4.781833  790.5614          [PC(32:1)+OAc][PE(35:1)+OAc][PS(36:0)+-H]
# 154  2.105127e-02 9.052746e-02  2.368832125 0.11643320 0.6601969 -5.105531  824.7721                                     [TG(48:0)+NH4]
# 80   2.236262e-01 3.825338e-01  2.260270329 0.12703673 0.6601969 -5.215086  762.6018                           [PC(34:0)+H][PE(37:0)+H]
# 213 -6.047636e-01 3.221100e-01 -2.210692987 0.13228967 0.6601969 -5.265777  950.8172                                     [TG(58:7)+NH4]
# 215 -8.504944e-02 4.865368e-02 -2.185338717 0.13508301 0.6601969 -5.291855  963.7433                                     [TG(59:10)+Na]
# 200 -2.297696e-01 1.240655e-01 -2.183406248 0.13529897 0.6601969 -5.293847  931.7710                          [TG(56:5)+Na][TG(58:8)+H]
# 201 -1.101939e+00 5.996685e-01 -2.155122797 0.13851052 0.6601969 -5.323068  926.8153                                     [TG(56:5)+NH4]
# 205 -5.860829e-01 3.231625e-01 -2.154031284 0.13863638 0.6601969 -5.324198  922.7853                                     [TG(56:7)+NH4]
# 189 -3.319595e-01 1.847236e-01 -2.153284148 0.13872262 0.6601969 -5.324972  905.7567                          [TG(54:4)+Na][TG(56:7)+H]
# 203 -1.987214e-01 1.089188e-01 -2.152383583 0.13882666 0.6601969 -5.325905  929.7560                          [TG(56:6)+Na][TG(58:9)+H]
# 192 -1.295042e+00 7.296724e-01 -2.145353384 0.13964223 0.6601969 -5.333191  898.7848                                     [TG(54:5)+NH4]





#########################################################################
################ Compare activated vs naive, in vitro ###################
#########################################################################

mscond_activated <- mscond
mscond_activated$tissue[mscond$tissue %in% c("th0","th1","th2")] <- "activated"

keep <- mscond_activated$tissue %in% c("activated","naive")
form <- ~1+mousenum+tissue
msde_activated <- doMStest(keep, form, mscond_orig = mscond_activated)
msde_activated

#sum(msde_activated$adj.P.Val<0.01)
msde_activated_denames <- msde_activated$name[msde_activated$adj.P.Val<0.05]
length(msde_activated_denames)


#########################################################################
################ Compare activated subsets ##############################
#########################################################################

keep <- mscond$tissue %in% c("th2","th0")
form <- ~1+mousenum+tissue
doMStest(keep, form)[1:20,]

keep <- mscond$tissue %in% c("th2","th0")
form <- ~1+tissue
doMStest(keep, form)[1:20,]


keep <- mscond$tissue %in% c("th1","th0")
form <- ~1+tissue
doMStest(keep, form)[1:20,]

keep <- mscond$tissue %in% c("th2","th1")
form <- ~1+tissue
doMStest(keep, form)[1:20,]


#########################################################################
################ Compare activated vs naive, in vitro ###################
#########################################################################


keep <- mscond$tissue %in% c("spleen","lung")
form <- ~1+tissue+isko
doMStest(keep, form)
# logFC      AveExpr            t     P.Value adj.P.Val         B       m.z                                               name
# 7   -7.497406e-02 4.035448e-01 -3.237509289 0.009436898 0.9991549 -3.693240  792.5769                       [PC(32:0)+OAc][PE(35:0)+OAc]
# 100 -5.583760e-03 1.103166e-02 -2.662399052 0.024688587 0.9991549 -4.667792  816.6491                           [PC(38:1)+H][PE(41:1)+H]
# 86  -4.854576e-03 1.339092e-02 -2.409277961 0.037795514 0.9991549 -5.089655  774.6012                           [PC(35:1)+H][PE(38:1)+H]
# 71  -7.671153e-02 4.745558e-01 -2.347143161 0.041949255 0.9991549 -5.191785  692.5230                           [PC(29:0)+H][PE(32:0)+H]
# 95  -3.273026e-02 1.867760e-01 -2.176893041 0.055753322 0.9991549 -5.467705  786.6027                           [PC(36:2)+H][PE(39:2)+H]
# 6    2.448940e-02 4.469553e-02  2.134997993 0.059773895 0.9991549 -5.534579  701.3966                                      [PA(34:5)+Cl]
# 8   -1.885053e-02 9.085114e-02 -2.091071514 0.064288416 0.9991549 -5.604201  790.5614          [PC(32:1)+OAc][PE(35:1)+OAc][PS(36:0)+-H]
# 21  -2.037345e-02 7.905048e-02 -2.020752089 0.072201482 0.9991549 -5.714519  738.5098                                      [PE(36:4)+-H]
# 10  -2.190185e-02 1.297840e-01 -1.795628473 0.104176946 0.9991549 -6.056726  816.5798                        [PC(34:2)+OAc][PS(38:1)+-H]

keep <- mscond$tissue %in% c("spleen","lung")
form <- ~1+isko
doMStest(keep, form)  ######## it appears that adding a factor for spleen and lung makes little difference

# logFC      AveExpr             t     P.Value adj.P.Val         B       m.z                                               name
# 7   -7.497406e-02 4.035448e-01 -3.1990267092 0.008897328 0.9992309 -3.666621  792.5769                       [PC(32:0)+OAc][PE(35:0)+OAc]
# 71  -7.671153e-02 4.745558e-01 -2.4639082193 0.032290018 0.9992309 -4.951808  692.5230                           [PC(29:0)+H][PE(32:0)+H]
# 6    2.448940e-02 4.469553e-02  2.2432614284 0.047395204 0.9992309 -5.323039  701.3966                                      [PA(34:5)+Cl]
# 100 -5.583760e-03 1.103166e-02 -2.1186253888 0.058733899 0.9992309 -5.527379  816.6491                           [PC(38:1)+H][PE(41:1)+H]
# 10  -2.190185e-02 1.297840e-01 -1.8867032186 0.086982290 0.9992309 -5.894315  816.5798                        [PC(34:2)+OAc][PS(38:1)+-H]
# 95  -3.273026e-02 1.867760e-01 -1.8312092550 0.095403444 0.9992309 -5.979100  786.6027                           [PC(36:2)+H][PE(39:2)+H]
# 80  -1.517889e-02 1.113383e-01 -1.7596946047 0.107356869 0.9992309 -6.086409  762.6018                           [PC(34:0)+H][PE(37:0)+H]

