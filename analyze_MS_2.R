#Acetoxy group, abbreviated AcO or OAc

# CE = Cholesteryl Esters
# LPC= Lysophosphocholine
# PA= phosphatic acid
# PS= phosphatidylglycerol   (PG on wikipedia)
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
# msdataposnames <- msdatapos[1,-1]
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





#########################################################################
############## Linear model to compare states  ##########################
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

############ The most DE in activation
# [1] "[TG(48:3)+NH4]"                         "[TG(49:0)+NH4]"                         "[TG(47:0)+NH4]"                        
# [4] "[LysoPC(18:1)+Na][LysoPC(20:4)+H]"      "[PC(40:5)+Na][PC(42:8)+H][PE(43:5)+Na]" "[PI(32:2)+Cl]"                         
# [7] "[PA(34:5)+Cl]"                          "[PI(32:2)+Cl]"     

#########################################################################
################ Compare activated subsets ##############################
#########################################################################

keep <- mscond$tissue %in% c("th2","th0")
form <- ~1+mousenum+tissue
v<-doMStest(keep, form)
v[v$name %in% msde_activated_denames,]
#[1:20,]

keep <- mscond$tissue %in% c("th2","th0")
form <- ~1+tissue
v<-doMStest(keep, form)
v[v$name %in% msde_activated_denames,]  #there is one DE 0.06

keep <- mscond$tissue %in% c("th1","th0")
form <- ~1+tissue
v<-doMStest(keep, form)
v[v$name %in% msde_activated_denames,]   #two have p<0.05


keep <- mscond$tissue %in% c("th2","th1")
form <- ~1+tissue
v<-doMStest(keep, form)
v[v$name %in% msde_activated_denames,]  #### 3 genes 0.01
### [PA(34:5)+Cl]     [PC(40:5)+Na][PC(42:8)+H][PE(43:5)+Na]     [PI(32:2)+Cl]


#########################################################################
################ Compare Rora KO ########################################
#########################################################################

keep <- mscond$tissue %in% c("spleen","lung")
form <- ~1+tissue+isko
v<-doMStest(keep, form)
v
v[v$name %in% msde_activated_denames,]
### [PA(34:5)+Cl]   is the only activation related one


#             logFC      AveExpr            t     P.Value adj.P.Val         B       m.z                                               name
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

#             logFC      AveExpr             t     P.Value adj.P.Val         B       m.z                                               name
# 7   -7.497406e-02 4.035448e-01 -3.1990267092 0.008897328 0.9992309 -3.666621  792.5769                       [PC(32:0)+OAc][PE(35:0)+OAc]
# 71  -7.671153e-02 4.745558e-01 -2.4639082193 0.032290018 0.9992309 -4.951808  692.5230                           [PC(29:0)+H][PE(32:0)+H]
# 6    2.448940e-02 4.469553e-02  2.2432614284 0.047395204 0.9992309 -5.323039  701.3966                                      [PA(34:5)+Cl]
# 100 -5.583760e-03 1.103166e-02 -2.1186253888 0.058733899 0.9992309 -5.527379  816.6491                           [PC(38:1)+H][PE(41:1)+H]
# 10  -2.190185e-02 1.297840e-01 -1.8867032186 0.086982290 0.9992309 -5.894315  816.5798                        [PC(34:2)+OAc][PS(38:1)+-H]
# 95  -3.273026e-02 1.867760e-01 -1.8312092550 0.095403444 0.9992309 -5.979100  786.6027                           [PC(36:2)+H][PE(39:2)+H]
# 80  -1.517889e-02 1.113383e-01 -1.7596946047 0.107356869 0.9992309 -6.086409  762.6018                           [PC(34:0)+H][PE(37:0)+H]


#### can say: all the major differences in Rora are of PC / PE type
# PC= phosphocholine
# PE= Phospoethanolamine
# TG/DG are glycerol. 
