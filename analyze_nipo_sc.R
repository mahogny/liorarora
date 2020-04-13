#######################################################################
################## Install packages ###################################
#######################################################################
if(FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  biocLite("topGO")
  biocLite("GO.db")
  biocLite("DESeq")
}


#set.seed(1234567890)
library(gplots)
library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(limma)
library(sqldf)
library(DESeq2)
library(DESeq)
library(xlsx)

source("common_functions.R")

#########
## This function hacks around problems in seurat. both me and zhichau has issues with the subsetting
reinit_seurat<-function(seu,sele){
  seu_new<-CreateSeuratObject(raw.data = seu@raw.data[,colnames(seu@raw.data)[sele]], 
                              min.cells = 2, project = "T cells", names.field = 2, names.delim = "\\-")
  for(feat in colnames(seu@meta.data)){
    seu_new <- AddMetaData(seu_new, setNames(seu@meta.data[sele,][[feat]], rownames(seu@meta.data[sele,])), feat)
  }
  seu_new
}

reinit_seurat_full<-function(seu_orig,sele){
  seu <- reinit_seurat(seu_orig, sele)
  #can't get below from the first object?
  seu <- NormalizeData(object = seu)
  seu@var.genes <- seu_orig@var.genes
  #seu <- FindVariableGenes(object = seu, x.low.cutoff = 0,x.high.cutoff = 5, y.cutoff = 0.1, do.plot = FALSE)
  seu <- ScaleData(object = seu, genes.use = seu@var.genes, model.use = "negbinom")
  seu
}


######################################################################
### Load data and perform QC #########################################
######################################################################

######## Load counts and metadata after QC
df <- readRDS('nipo_ss2/df1.rds')
phn <- df$phn
# QC already done

seu<-CreateSeuratObject(raw.data = df$cts, min.cells = 2, project = "T cell", 
                        names.field = 2, names.delim = "\\-")
for(feat in colnames(df$phn)){
  seu <- AddMetaData(seu, setNames(df$phn[[feat]], rownames(df$phn)), feat)
}

##### Normalize and find variable genes
seu <- NormalizeData(object = seu)
seu <- FindVariableGenes(object = seu, x.low.cutoff = 0,x.high.cutoff = 5, y.cutoff = 0.1)
length(x = seu@var.genes)

##### Add TraCeR annotation: NK-cells
df <- read.csv("nipo_ss2/NK.txt",header = F)
nk <- rep(0,dim(seu@raw.data)[2])
nk[colnames(seu@raw.data) %in% df$V1 ]<-1
seu <- AddMetaData(seu, setNames(nk, colnames(seu@raw.data)), "NK")

##### Add TraCeR annotation: doublets in library
df <- read.csv("nipo_ss2/doublet-like.txt",header = F)
dl <- rep(0,dim(seu@raw.data)[2])
dl[colnames(seu@raw.data) %in% df$V1 ]<-1
seu <- AddMetaData(seu, setNames(dl, colnames(seu@raw.data)), "doublet")


######################################################################
### Check exon4 usage ################################################
######################################################################
#The exon4 (ENSMUSE00001263945) was knock-out by loxp system. 
#It is a conditional KO, so the KO is effective only when the cell DO NOT express cd4. 

#GGCTTTTTCAGGAGAAGTCAGCAGAGCAATGCCACCTACTCCTGTCCTCGTCAGAAGAACTGTTTGATTGATCGGACCAGCAGAAACCGCTGCCAGCATTGTCGGCTGCAGAAATGCCTGGCCGTGGGGATGTCTCGAGATG


df_ex <- read.csv('nipo_ss2/exon35.csv', row.names = 1)
colnames(df_ex)<-gsub('\\.','#',gsub('X','',colnames(df_ex)))
seu <- AddMetaData(seu, setNames(as.integer(df_ex["ENSMUSE00001263945",]), colnames(df_ex)), "exon4")
exon4_rel <- log(1+seu@meta.data$exon4/(1+seu@raw.data["Rora",]))
seu <- AddMetaData(seu, setNames(exon4_rel, colnames(seu@raw.data)), "exon4_rel")
seu <- AddMetaData(seu, setNames(seu@meta.data$exon4>0,   colnames(seu@raw.data)), "exon4_pos")

###### Plot relative count
p<-ggplot(seu@meta.data, aes(genotype, exon4_rel))+
  geom_jitter(aes(colour = genotype), height = 0, width = 0.3)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")
plot(p)
ggsave("out.niposs2/ko_exon4_relcount.pdf",p)

###### Plot raw count
p<-ggplot(seu@meta.data, aes(genotype, exon4))+
  geom_jitter(aes(colour = genotype), height = 0, width = 0.3)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")
plot(p)
ggsave("out.niposs2/ko_exon4_rawcount.pdf",p)


######################################################################
### Annotate particular expression levels ############################
######################################################################
seu <- AddMetaData(seu, setNames(seu@raw.data["Rora",],   colnames(seu@raw.data)), "Rora_exprs")
seu <- AddMetaData(seu, setNames(seu@raw.data["Rora",]>1, colnames(seu@raw.data)), "Rora_exprs_pos")
#seu <- AddMetaData(seu, setNames(seu@raw.data["Rora",]>5, colnames(seu@raw.data)), "Rora_exprs_pos")

#hist(seu@raw.data["Rora",],breaks=50)
#hist(seu@scale.data["Rora",],breaks=50)
#sum(seu@meta.data$Rora_exprs_pos)

seu <- AddMetaData(seu, setNames(seu@raw.data["Cd4",]>1, colnames(seu@raw.data)),  "Cd4_exprs_pos")
seu <- AddMetaData(seu, setNames(seu@raw.data["Sell",]>1, colnames(seu@raw.data)), "Sell_exprs_pos")

seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="lungs",  colnames(seu@raw.data)), "is_lung")
seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="spleen", colnames(seu@raw.data)), "is_spleen")
seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="MedLN",  colnames(seu@raw.data)), "is_MedLN")
seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="MLN",    colnames(seu@raw.data)), "is_MesLN")




######################################################################
### Cluster cells together ###########################################
######################################################################
set.seed(1234567890)
seu <- ScaleData(object = seu, genes.use = seu@var.genes, model.use = "negbinom")
seu <- RunPCA(object = seu, pc.genes = seu@var.genes, pcs.compute = 50, pcs.print = 1:2, maxit = 500, weight.by.var = FALSE)
seu <- RunTSNE(object = seu, dims.use = 1:25, do.fast = TRUE)
seu <- FindClusters(object = seu, reduction.type = "pca",  dims.use = 1:25, resolution = 0.5, save.SNN = TRUE, force.recalc = T )
clusterident <- seu@ident

###### Assign names
clusterident_mod <- clusterident
levels(clusterident_mod)[0+1] <- "Th2_day_30" #ok
levels(clusterident_mod)[1+1] <- "Naive_day_7_30" #ok
levels(clusterident_mod)[2+1] <- "Treg_day_7_30"
levels(clusterident_mod)[3+1] <- "Th1_day_30"  #changed
levels(clusterident_mod)[4+1] <- "Naive_day_30"  #changed
levels(clusterident_mod)[5+1] <- "Th1_day_7"
levels(clusterident_mod)[6+1] <- "noidea_day_7"
levels(clusterident_mod)[7+1] <- "Tfh"
#levels(clusterident_mod)[8+1] <- "Macrophages"  #disappeared as a group

#Hack in NK cells as a subset
clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "NK"))
clusterident_mod[seu@meta.data$NK==1] <- "NK"
seu@ident <- clusterident_mod

###### For easy referencing
subsets_cluster <- c("Th1_day_30","Th2_day_30","Treg_day_7_30","Naive_day_7_30")
subsets_name    <- c("Th1",       "Th2",       "Treg",         "Naive")


###### Plot the clustering
pdf("out.niposs2/clustering_allcell.pdf")
TSNEPlot(object = seu, do.label = TRUE)
dev.off()


###### Plot the clustering: Only CD4, no doublets
pdf("out.niposs2/clustering_cd4.pdf")
selcell <- colnames(seu@raw.data)[seu@meta.data$Cd4_exprs_pos & seu@meta.data$doublet<1]
seu_tsne <- RunTSNE(object = seu, dims.use = 1:25, do.fast = TRUE, cells.use = selcell)
TSNEPlot(object = seu_tsne, do.label = TRUE, cells.use = selcell)
dev.off()


###### Plot the clustering: Color by day
pdf("out.niposs2/clustering_allcell_byday.pdf")
thecol <- rep("blue",1670)
thecol[seu@meta.data$days.after.infection==7] <- "red"
thecol[seu@meta.data$days.after.infection==30] <- "green"
TSNEPlot(object = seu, colors.use=thecol)
dev.off()

# 
# ###### Various genes
# pdf("out.niposs2/clustering_allcell_markers.pdf",h=20,w=20)
# FeaturePlot(seu, c(
#   "Tbx21", "Gata3","Il4","Il13","Ifng", #Th1-2 axis
#   "Foxp3",       #Tregs
#   "Sell","Cd44", #Activated T cells
#   "Rorc","Rora",
#   "Il10","Cxcr5","Ccr3",
#   "Cd4",
#   #"Il1b",
#   #"Ear2",
#   
#   "Cd19",
#   #Cd19, Cd20, Cd34, Cd38, and Cd45r  #immature B cell
#   #"Itgax",
#   #"Col14a1","Trim47","Ccl12","Ccl7","Apoe","Mgl2","Cd63","Ly86","Cxcl16","Ccl8","Ccl2","Anxa3","Ccl3",
#   "Cd68",                         #pulls out macrophages
#   #"C1qa","C1qb","C1qc","Sirpa",  #pulls out macrophages
#   #"Il17a","Il17f","Il21","Il22",
#   #"Ptprc", #Cd45
#   #"Cd8a",
#   #"Cd14","Lyz2","Ms4a7",
#   
#   #http://www.abcam.com/primary-antibodies/effector-t-cell-markers
#   
#   # 	CD8+, CD45RA+, CD45RO-, CCR7+, CD28+	IFN-γ+, IL-2+    Naive Cd8
#   #   CD8+, CCR7-	IFN-γ+, Perforin+, Granzyme+    cytotoxic cd8
#   # 	CD4+, CCR3+, CCR6+	IL-9+
#   #	CD4+, CCR6+, CCR4+, NK1.1+	IL-17+
#   #	CD4+, CCR10+, CCR4+, CCR6+ Il22   Th22
#   "Il23r",#gammadelta
# #  "Il21","Bcl6","Cxcr5","Cd40lg","Icos",  #Tfh markers
# 
#   "Il1rl1",   #nuocyte markers, these are mixed with the Th2   IL-17BR, ICOS or ST2=il1rl1
#   #IL-17A, IL-17F, IL-21, and IL-22
#   "days.after.infection", "is_lung","is_spleen","is_MedLN","is_MesLN"), 
#   cols.use = c("green", "blue"),pt.size = 1)
# dev.off()
# 
# 
# pdf("out.niposs2/clustering_allcell_markers_2.pdf",h=20,w=20)
# selcell <- colnames(seu@raw.data)[seu@meta.data$Cd4_exprs_pos & seu@meta.data$doublet<1 & #seu@meta.data$genotype!="control mouse" &
#                                     seu@meta.data$days.after.infection>=30]
# FeaturePlot(seu, c(
#   #"Fhl2","Gzmc","Tnfsf14","Cd8b1","Map9","Gpr55","Il10",
#   "Gzmc","Fhl2"
#   # "Cndp2","Gnpda1","Mad2l1bp","Cxcr5","Bpgm","Il10"
#   # ,"Rora","Ppara","Pparg","Nr1d1","Nr1d2","Rxra","Vdr"
#   # ,"Rorb","Rorc"
#   ), 
#   cols.use = c("green", "blue"),pt.size = 2, cells.use = selcell)
# dev.off()
# #hist(seu@scale.data["Cxcr5",],breaks=100)



# Effector   B cell?
#grep("Cd45",ensconvert$mgi_symbol)


################ Markers as violins ####################
pdf("out.niposs2/clustering_allcell_violin.pdf",h=20,w=20)
VlnPlot(seu, c("Gata3","Sell","Tbx21","Cd4","Foxp3","Il10"),x.lab.rot = TRUE)
dev.off()



pdf("out.niposs2/markers_rora_il10.pdf",h=5,w=9)
VlnPlot(seu, c("Rora","Il10"),x.lab.rot = TRUE, point.size.use = 0.5)
dev.off()

pdf("out.niposs2/markers_liorafav.pdf",h=30,w=15)
VlnPlot(seu,
        c(
          "Tnfrsf25",
          "Il2",
          "S100a4",
          "Actn2",
          "Arntl",
          "Nebl",
          "Cxcr6",
          "Fasl",
          "Tbx21",
          "Gzma",
          "Hhex",
          "Stx2",
          "Aldh7a1") 
        ,x.lab.rot = TRUE, point.size.use = 0.5)
dev.off()






selcell <- colnames(seu@raw.data)[seu@meta.data$Cd4_exprs_pos & seu@meta.data$doublet<1 & 
                                    seu@meta.data$days.after.infection>=30 & clusterident_mod==idtreg]
VlnPlot(seu, c("Il10"),x.lab.rot = TRUE, cells.use=selcell)


######################################################################
### Subset CD4+ cells based on known markers #########################
######################################################################

### Get marker genes from Th-express. 
d <-         read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_naive.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE]
d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_Th1.txt",  sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_Th17.txt", sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_iTreg.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_nTreg.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
thep <- d[,c(2,4,6,8,10)]
thefc <- d[,c(1,3,5,7,9)]
colnames(thep)<-c("Naive/Th0","Th1","Th17","iTreg","nTreg")
colnames(thefc)<-c("Naive/Th0","Th1","Th17","iTreg","nTreg")

for(i in 1:ncol(thep)){
  thep[is.na(thep[,i]),i] <- 1
}
for(i in 1:ncol(thefc)){
  thefc[is.na(thefc[,i]),i] <- 0
}
thmarkergenes <- names(which(apply(thep<1e-17 & abs(thefc)>5,1,any)))
thmarkergenes_sym <- sort(togenesym(thmarkergenes))
#thmarkergenes
# c("Foxp3","Tbx21","Gata3") %in% thmarkergenes_sym
# length(thmarkergenes_sym)



######################################################################
### Compare: infected KO CD4+ vs infected control CD4+, ##############
### for each cell type (Th1, Th2, Treg)                 ##############
######################################################################


#> CD4+ clusters, treatment: infected, exp: KO, days.after.infection:30, 
#  Rora>1, exon4 in KO <1, tissue: lung

getdeone_sc <- function(theid="", onlyrora=FALSE){
  set.seed(123)
  ### Split cells into 2 groups (and neither)
  newident <- rep(0,length(seu@ident))
  newident[!seu@meta.data$exon4_pos &
             seu@meta.data$Cd4_exprs_pos & 
             seu@meta.data$days.after.infection>=30 &
             !seu@meta.data$exon4_pos &
             seu@meta.data$genotype=="rora (CD4) KO"] <- "ko"
  newident[seu@meta.data$Cd4_exprs_pos & 
             seu@meta.data$days.after.infection>=30 &
             seu@meta.data$genotype=="control mouse"] <- "ctrl"
  if(onlyrora){
    newident[!seu@meta.data$Rora_exprs_pos] <- 0
  }
  if(theid!=""){
    newident[clusterident_mod!=theid] <- 0
  }
  newident <- as.factor(newident)
  print(sum(newident=="ko"))
  print(sum(newident=="ctrl"))
  names(newident) <- names(seu@ident)
  seu@ident <- as.factor(newident)
  
  ### Perform the DE. for this number ideally with DEseq
  #guess I need to make an object with rounded read counts?
  one_de <- FindMarkers(seu, ident.1 = "ko", ident.2 = "ctrl")#, test.use = "DESeq2")
  #one_de <- FindMarkers(seu, ident.1 = "ko", ident.2 = "ctrl", test.use = "DESeq2")
  one_de  
}

subsets_cluster <- c("Th1_day_30","Th2_day_30","Treg_day_7_30","Naive_day_7_30")
subsets_name <- c("Th1","Th2","Treg","Naive")

seu@raw.data <- round(seu@raw.data)

#de_sc_deseq <- de_sc

de_sc <- list()
#de_sc[["all"]] <- getdeone_sc("")
for(i in 1:length(subsets_name)){
  print(subsets_cluster[i])
  de_sc[[subsets_name[i]]] <- getdeone_sc(subsets_cluster[i])
}

# de_sc_rora <- list()
# de_sc_rora[["all"]] <- getdeone_sc("", onlyrora=TRUE)
# for(i in 1:length(subsets_name)){
#   print(subsets_cluster[i])
#   de_sc_rora[[subsets_name[i]]] <- getdeone_sc(subsets_cluster[i], onlyrora=TRUE)
# }




