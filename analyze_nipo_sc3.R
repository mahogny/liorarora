
## redone analysis for paper revision


#######################################################################
################## Install packages ###################################
#######################################################################
if(FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  biocLite("topGO")
  biocLite("GO.db")
  biocLite("DESeq")
  #install.packages('devtools')
  # Replace '2.3.0' with your desired version
  devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))
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

# seu<-CreateSeuratObject(counts = df$cts, min.cells = 2, project = "T cell",
#                         names.field = 2, names.delim = "\\-")
seu<-CreateSeuratObject(raw.data = df$cts, min.cells = 2, project = "T cell",
                         names.field = 2, names.delim = "\\-")
for(feat in colnames(df$phn)){
  seu <- AddMetaData(seu, setNames(df$phn[[feat]], rownames(df$phn)), feat)
}

##### Normalize and find variable genes
seu <- NormalizeData(object = seu)
seu <- FindVariableGenes(object = seu, x.low.cutoff = 0,x.high.cutoff = 5, y.cutoff = 0.1)
#seu <- FindVariableFeatures(object = seu)
#length(x = seu@var.genes)

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


########################################
### Build table for arrayexpress
dat_sc_plate <- read.csv("plate_sc/liora_pools_red.csv")
xt_sub <- read.csv("plate_sc/xt_map.csv")
dat_sc_plate <- merge(dat_sc_plate, xt_sub)
dat_sc_plate$file1 <- sprintf("22719_%s_%s.1.fq.gz",dat_sc_plate$lane, dat_sc_plate$library_num)
dat_sc_plate$file2 <- sprintf("22719_%s_%s.2.fq.gz",dat_sc_plate$lane, dat_sc_plate$library_num)
dat_sc_plate$rname <- sprintf("22719_%s#%s",        dat_sc_plate$lane, dat_sc_plate$library_num)
dat_sc_plate$pname <- sprintf("Th nipo_%s_%s",        dat_sc_plate$lane, dat_sc_plate$library_num)
dat_sc_plate <- dat_sc_plate[dat_sc_plate$library_num<=380,]
rownames(dat_sc_plate) <- dat_sc_plate$rname
# sum(dat_sc_plate$rname %in% rownames(seu@meta.data))
# dim(dat_sc_plate)  
#write.csv(dat_sc_plate, "plate_sc/for_arrayexpress.csv")



tofix <- read.csv("plate_sc/sdrf.tsv", sep="\t",stringsAsFactors = TRUE, as.is = TRUE, check.names = FALSE)
tofix$`Comment[LIBRARY_STRATEGY]` <-"RNA-Seq"
tofix$`Comment[NOMINAL_LENGTH]` <-"150"
tofix$`Comment[LIBRARY_STRAND]` <-"not applicable"
tofix$`Comment[NOMINAL_SDEV]` <-30
tofix$`Comment[LIBRARY_SOURCE]` <-"TRANSCRIPTOMIC SINGLE CELL"
tofix$`Comment[LIBRARY_SELECTION]` <-"Oligo-dT"
tofix$`Comment[LIBRARY_LAYOUT]` <- "PAIRED"
tofix$`Array Data File` <- tofix$`Assay Name`
tofix$`Array Data File` <- str_replace_all(tofix$`Array Data File`,"Th nipo_","22719_")  
tofix1 <- tofix
tofix2 <- tofix
tofix1$`Array Data File` <- sprintf("%s.1.fastq.gz",tofix1$`Array Data File`)
tofix2$`Array Data File` <- sprintf("%s.2.fastq.gz",tofix2$`Array Data File`)


tofix_tot <- rbind(tofix1,tofix2)

write.table(tofix_tot, "plate_sc/sdrf.new.tsv")

######################################################################
### Check exon4 usage ################################################
######################################################################
#The exon4 (ENSMUSE00001263945) was knock-out by loxp system. 
#It is a conditional KO, so the KO is effective only when the cell DO NOT express cd4. 

#GGCTTTTTCAGGAGAAGTCAGCAGAGCAATGCCACCTACTCCTGTCCTCGTCAGAAGAACTGTTTGATTGATCGGACCAGCAGAAACCGCTGCCAGCATTGTCGGCTGCAGAAATGCCTGGCCGTGGGGATGTCTCGAGATG


df_ex <- read.csv('nipo_ss2/exon35.csv', row.names = 1)
colnames(df_ex)<-gsub('\\.','#',gsub('X','',colnames(df_ex)))


apply_exon_count <- function(seu){
  seu <- AddMetaData(seu, setNames(as.integer(df_ex["ENSMUSE00001263945",]), colnames(df_ex)), "exon4")
  exon4_rel <- log(1+seu@meta.data$exon4/(1+seu@raw.data["Rora",]))
  seu <- AddMetaData(seu, setNames(exon4_rel, colnames(seu@raw.data)), "exon4_rel")
  seu <- AddMetaData(seu, setNames(seu@meta.data$exon4>0,   colnames(seu@raw.data)), "exon4_pos")
}

seu <- apply_exon_count(seu)


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

addAnno <- function(seu){
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
}
seu <- addAnno(seu)



######################################################################
### Cluster all cells together #######################################
######################################################################
set.seed(1234567890)
seu <- ScaleData(object = seu, genes.use = seu@var.genes, model.use = "negbinom")
seu <- RunPCA(object = seu, pc.genes = seu@var.genes, pcs.compute = 50, pcs.print = 1:2, maxit = 500, weight.by.var = FALSE)
seu <- RunTSNE(object = seu, dims.use = 1:25, do.fast = TRUE)
seu <- FindClusters(object = seu, reduction.type = "pca",  dims.use = 1:25, resolution = 0.5, save.SNN = TRUE, force.recalc = T )
clusterident <- seu@ident
clusterident_nomarker <- seu@ident

### Original cluster names
clusterident_mod <- clusterident
seu@ident <- clusterident_mod
TSNEPlot(object = seu, do.label = TRUE)



##### identiy clusters
VlnPlot(seu,
        c(
          "Cd4",
          "Ptprc",   #"Cd45a",
#          "Itgal",   #Cd11a
 #         "Itgam",   #Cd11b
#          "Itgax",   #"Cd11c",  #cluster 6
          "Itgae",  #Cd103    mainly in treg
          "Rora",
          "Sell",

          ### DC markers
          #"Flt3",
          "Btla",   #cluster 7, and b cell
          #"Ly6c",
          #"Cd1c"
          #"Bst2",
          #"Tlr7"

          ### Th1 marker
          #"Tbx21"


          ### Treg marker
          #"Foxp3",  

          ### Th17 marker
          #"Rorc"
          #"Il17a",  #not there

          ### Th2 marker
          #"Gata3",
          #"Il4"

          #"Cxcr6"  #Th17 marker??  0,4,5,7. especially 5
          #"Cd8a",   #nowhere, or spread out
          #"Cd66b", #"Cd66b",  #granulocyte  -- removed, no exp?
#          "Itgam"  #Cd11b

          #### B cell markers
          #"Ms4a1"   #B cell marker  "Cd20"   #big in cluster 4   
          #"Cd19",   #big in cluster 4

          #### Macrophage markers: CD14, CD16, CD64, CD68, CD71 and CCR5;
          #"Cd14",   #macrophage marker   in cluster 6 if any
          #"Fcgr3a","Cd16"
          #"Fcgr1"#,"Cd64*"     ## cluster 6
        ) 
        ,x.lab.rot = TRUE, point.size.use = 0.5)
### with 0-7, pcs 50, 1:25
#3=th1, tbx21 high  
#0=th2 (gata3)   mainly d30
#2=Treg (foxp3)  all days
#4=day 7 (sell+)     B cells
#5=NK
#6   weird sell+  day 7 (small)
#7   weird sell+  day 7 (small)    ptprc+
#1   also naive?  sell+ ptprc+ 

TSNEPlot(object = seu, do.label = TRUE)
#6,7,4 weird  -- days?
#5 is NK cells, based on annotation below
seu@meta.data$days.after.infection



#Hack in NK cells as a subset
clusterident_mod <- clusterident
clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "NK"))
clusterident_mod[seu@meta.data$NK==1] <- "NK"
seu@ident <- clusterident_mod

#Plot by day
colorby_day <- function(seu, clusterident) {
  clusterident_mod <- clusterident
  clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "d30"))
  clusterident_mod[seu@meta.data$days.after.infection==30] <- "d30"
  clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "d0"))
  clusterident_mod[seu@meta.data$days.after.infection==0] <- "d0"
  clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "d7"))
  clusterident_mod[seu@meta.data$days.after.infection==7] <- "d7"
  seu@ident <- clusterident_mod
  seu
}
#TSNEPlot(object = seu, do.label = TRUE)

######### Need to save a picture of just this?

######## Color by genotype
colorby_genotype <- function(seu, clusterident) {
  clusterident_mod <- clusterident
  clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "HETS"))
  clusterident_mod[seu@meta.data$genotype=="RoraTEAL B6 HETS"] <- "HETS"
  clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "KO"))
  clusterident_mod[seu@meta.data$genotype=="rora (CD4) KO"] <- "KO"
  clusterident_mod <- factor(clusterident_mod, levels = c(levels(clusterident_mod), "control"))
  clusterident_mod[seu@meta.data$genotype=="control mouse"] <- "control"
  seu@ident <- clusterident_mod
  seu
}
#TSNEPlot(object = seu, do.label = TRUE)



###### Assign names
clusterident_mod <- clusterident
levels(clusterident_mod)[0+1] <- "Th2"
levels(clusterident_mod)[1+1] <- "Naive"
levels(clusterident_mod)[2+1] <- "Treg"
levels(clusterident_mod)[3+1] <- "Th1"
levels(clusterident_mod)[4+1] <- "B cell"
levels(clusterident_mod)[5+1] <- "NK"
levels(clusterident_mod)[6+1] <- "Macrophage"
levels(clusterident_mod)[7+1] <- "DCs"
seu@ident <- clusterident_mod
seu@meta.data$inferred_celltype <- clusterident_mod
TSNEPlot(object = seu, do.label = TRUE)


###### For easy referencing
subsets_cluster <- c("Th1","Th2","Treg","Naive")
subsets_name    <- c("Th1","Th2","Treg","Naive")

TSNEPlot(object = seu, do.label = TRUE)


######################################################################
### Cluster only T helper cells together #############################       
######################################################################

#Only keep Th cells
keep <- seu@meta.data$inferred_celltype %in% subsets_cluster
seu_sub <- SubsetData(seu, cells.use = keep)


set.seed(1234567890)
seu_sub <- ScaleData(object = seu_sub, genes.use = seu_sub@var.genes, model.use = "negbinom")
seu_sub <- RunPCA(object = seu_sub, pc.genes = seu_sub@var.genes, pcs.compute = 50, pcs.print = 1:2, maxit = 500, weight.by.var = FALSE)
seu_sub <- RunTSNE(object = seu_sub, dims.use = 1:25, do.fast = TRUE)
clusterident_sub <- seu_sub@ident

TSNEPlot(object = seu_sub, do.label = TRUE)   #By cell type
TSNEPlot(object = colorby_genotype(seu_sub, clusterident_sub), do.label = TRUE)
TSNEPlot(object = colorby_day(seu_sub, clusterident_sub), do.label = TRUE)



###### 
plot_colorby_tissue <- function(seu) {
  thecol <- rep("gray",nrow(seu@meta.data))
  thecol[as.character(seu@meta.data$tissue)=="lungs"]  <- "red"
  thecol[as.character(seu@meta.data$tissue)=="spleen"] <- "green"
  thecol[as.character(seu@meta.data$tissue)=="MedLN"]  <- "blue"
  thecol[as.character(seu@meta.data$tissue)=="MLN"]    <- "blue"  #Why both?
  plot(
    seu@dr$tsne@cell.embeddings[,1],
    seu@dr$tsne@cell.embeddings[,2],
    col=thecol, pch=19,cex=0.5,
    xlab="tSNE_1",ylab="tSNE_2",main="tissue"
  )
}
pdf("paper/revision fig/nipo tissue.pdf")
plot_colorby_tissue(seu_sub)
dev.off()

### celltype as color, rora as shape - impossible to read
plot_gene_exp_over_celltype <- function(seu, gene="Rora"){
  thecol <- rep("gray",nrow(seu@meta.data))
  thecol[as.character(seu@ident)=="Naive T cell"]  <- "red"
  thecol[as.character(seu@ident)=="Th1"] <- "green"
  thecol[as.character(seu@ident)=="Th2"]  <- "blue"
  thecol[as.character(seu@ident)=="Treg"]    <- "purple"  
  
  thepch <- rep(19, nrow(seu@data))
  #thepch[seu@data[gene,]>1] <- 20

  plot(
    seu@dr$tsne@cell.embeddings[,1],
    seu@dr$tsne@cell.embeddings[,2],
    col=thecol, pch=thepch,cex=0.5,
    xlab="tSNE_1",ylab="tSNE_2", main=gene
  )
}
pdf("paper/revision fig/nipo celltype.pdf")
plot_gene_exp_over_celltype(seu_sub)
dev.off()

plot_gene_exp <- function(seu, gene="Rora", binarize=F){
  thecol <- rgb(seu@data[gene,]/max(seu@data[gene,]),0,0)
  if(binarize){
    thecol <- rgb(round(seu@data[gene,]>1),0,0)
  }
    
  plot(
    seu@dr$tsne@cell.embeddings[,1],
    seu@dr$tsne@cell.embeddings[,2],
    col=thecol, pch=19,cex=0.5,
    xlab="tSNE_1",ylab="tSNE_2", main=gene
  )
}
plot_gene_exp(seu_sub,"Rora")
pdf("paper/revision fig/nipo rora binary.pdf")
plot_gene_exp(seu_sub,"Rora",T)
dev.off()
plot_gene_exp(seu_sub,"Tnfrsf25")
plot_gene_exp(seu_sub,"Gata3")
plot_gene_exp(seu_sub,"S100a4")
#plot_gene_exp(seu_sub,"Alox8")
plot_gene_exp(seu_sub,"Il1rl1")  #Il33 receptor
plot_gene_exp(seu_sub,"Tgfb1")
plot_gene_exp(seu_sub,"Il13")


plot_cyto_exp <- function(seu){
  norm_max <- function(g1) {
    if(max(g1)==0)
      g1
    else
      g1/max(g1)
  }
  #"Tgfb1"  behaves badly
  genes <- c("Il4","Ifng","Il17a","Il5","Il10","Il13","Il2")
  #g1 <- colSums(seu@data[genes,])
  g1 <- apply(seu@data[genes,],2,norm_max)
  print(g1)
  g1 <- colMaxs(g1)
#  g1 <- colSums(g1)

  thecol <- rgb(g1/max(g1),0,0)
  plot(
    seu@dr$tsne@cell.embeddings[,1],
    seu@dr$tsne@cell.embeddings[,2],
    col=thecol, pch=19,cex=0.5,
    xlab="tSNE_1",ylab="tSNE_2", main="Cytokines"
  )
}
plot_cyto_exp(seu_sub)



##### Violins of interesting genes
VlnPlot(seu_sub,
        c(
          #"Rora",
          #"Sell",
          
          #From activation section of paper
          #"Arntl",
          #"Tnfrsf25",
          # "Alox8",
          # "S100a4",
          # "Cxcr6",
          # "St6galnac3",
          # "Emb",
          
          
          ### Receptors
          # "Il1rl1",   #IL33 receptor    Th2/Treg
          # "Ccr5",     #activated
          # "Ccr3",     #not really
          # "Il15ra",   #all uniformly
          # "Il6ra",    #all but Th1   -- mainly naive
          # "Ccr4"      #Th2/Treg
          
          ### Final model, most relevant cytokines  IL33 (to Il1rl1), SDF1a (Cxcr4), CCL7 (many!), CCL22 (Ccr4), 
          # "Il1rl1",   #IL33 receptor    Th2/Treg
          # "Cxcr4",      #a bit everywhere
          # "Ccr4"      #Th2/Treg but also some in naive
             
          
        ) 
        ,x.lab.rot = TRUE, point.size.use = 0.5)



### To highlight, final list
pdf("paper/revision fig/highlight_violins.pdf", h=6)
VlnPlot(seu_sub,
        c(
          "Rora",
          
          "Ccr5",
          "Cxcr4",
          "Cxcr6",
          "Il1rl1",
          "Il6ra",
          
          "Arntl",
          "Alox8",
          "S100a4",
          "Tnfrsf25"

        ) 
        ,x.lab.rot = TRUE, point.size.use = .3,size.title.use = 8, size.x.use = 8)
dev.off()




######################################################################
### Cluster only T helper cells together #############################       OLD
######################################################################

#Only keep Th cells
keep <- seu@meta.data$ident %in% subsets_cluster

df_th <- readRDS('nipo_ss2/df1.rds')
phn <- df_th$phn
df_th$cts <- df_th$cts[,keep]
df_th$phn <- df_th$phn[keep,]
seu<-CreateSeuratObject(raw.data = df_th$cts, min.cells = 2, project = "T cell", 
                        names.field = 2, names.delim = "\\-")
for(feat in colnames(df_th$phn)){
  seu <- AddMetaData(seu, setNames(df_th$phn[[feat]], rownames(df_th$phn)), feat)
}

##### Normalize and find variable genes
seu <- NormalizeData(object = seu)
seu <- FindVariableGenes(object = seu, x.low.cutoff = 0,x.high.cutoff = 5, y.cutoff = 0.1)
length(x = seu@var.genes)

##### Recluster
set.seed(1234567891)
seu <- ScaleData(object = seu, genes.use = seu@var.genes, model.use = "negbinom")
seu <- RunPCA(object = seu, pc.genes = seu@var.genes, pcs.compute = 50, pcs.print = 1:2, maxit = 500, weight.by.var = FALSE)
seu <- RunTSNE(object = seu, dims.use = 1:25, do.fast = TRUE)
seu <- FindClusters(object = seu, reduction.type = "pca",  dims.use = 1:25, resolution = 0.5, save.SNN = TRUE, force.recalc = T )
seu <- addAnno(seu)
seu <- apply_exon_count(seu)
clusterident <- seu@ident

### Annotate these clusters
clusterident_mod <- clusterident
levels(clusterident_mod)[0+1] <- "Th2"
levels(clusterident_mod)[1+1] <- "Treg"
levels(clusterident_mod)[2+1] <- "Naive"
levels(clusterident_mod)[3+1] <- "Th1"  
levels(clusterident_mod)[4+1] <- "remove"
seu@ident <- clusterident_mod

# New "shortcuts"
subsets_cluster <- c("Th1","Th2","Treg","Naive")
subsets_name    <- c("Th1",       "Th2",       "Treg",         "Naive")

###### Plot the clustering. Exclude cluster 4. Focus on cluster names
pdf("out.niposs2/new_clustering_thcell.pdf")
TSNEPlot(object = seu, do.label = TRUE, cells.use = rownames(seu@meta.data)[clusterident_mod!="remove"])
dev.off()

###### Plot the clustering. Exclude cluster 4. Focus on tissue sources
thecol <- rep("gray",nrow(seu@meta.data))
thecol[as.character(seu@meta.data$tissue)=="lungs"]  <- "red"
thecol[as.character(seu@meta.data$tissue)=="spleen"] <- "green"
thecol[as.character(seu@meta.data$tissue)=="MedLN"]  <- "blue"
thecol[as.character(seu@meta.data$tissue)=="MLN"]    <- "blue"  #Why both?
TSNEPlot(object = seu, colors.use=thecol[clusterident_mod!="remove"],do.label = TRUE, cells.use = rownames(seu@meta.data)[clusterident_mod!="remove"])
pdf("out.niposs2/clustering_celltissue.pdf",width = 5.6, height = 7.6)
plot(
  seu@dr$tsne@cell.embeddings[,1][clusterident_mod!="remove"],
  seu@dr$tsne@cell.embeddings[,2][clusterident_mod!="remove"],
  col=thecol[clusterident_mod!="remove"], pch=19,cex=0.5,
  xlab="tSNE_1",ylab="tSNE_2"
)
dev.off()


#### Should plot Rora level!
thecol <- rgb(seu@data["Rora",]/max(seu@data["Rora",]),0,0)
pdf("out.niposs2/new_clustering_rora.pdf",width = 5.6, height = 7.6)
plot(
  seu@dr$tsne@cell.embeddings[,1][clusterident_mod!="remove"],
  seu@dr$tsne@cell.embeddings[,2][clusterident_mod!="remove"],
  col=thecol[clusterident_mod!="remove"], pch=19,cex=0.5,
  xlab="tSNE_1",ylab="tSNE_2"
)
dev.off()

# 
# ###### Plot the clustering: Color by day
# pdf("out.niposs2/clustering_allcell_byday.pdf")
# thecol <- rep("blue",nrow(seu@meta.data))
# thecol[seu@meta.data$days.after.infection==7] <- "red"
# thecol[seu@meta.data$days.after.infection==30] <- "green"
# plot(
#   seu@dr$tsne@cell.embeddings[,1],
#   seu@dr$tsne@cell.embeddings[,2],
#   col=thecol, pch=19
# )
# #TSNEPlot(object = seu, colors.use=thecol)
# dev.off()
# 


################ Markers as violins ####################
pdf("out.niposs2/clustering_thcell_violin.pdf",h=20,w=20)
VlnPlot(seu, c("Gata3","Sell","Tbx21","Cd4","Foxp3","Il10"),x.lab.rot = TRUE)
dev.off()



# pdf("out.niposs2/markers_rora_il10.pdf",h=5,w=9)
# VlnPlot(seu, c("Rora","Il10"),x.lab.rot = TRUE, point.size.use = 0.5)
# dev.off()

pdf("out.niposs2/violin_exgenes.pdf",h=8,w=15)
VlnPlot(seu,
        c(
          "Arntl",
          "Alox8",
          "Ccr4",
          "Cxcr6",
          "Fasl",
          "Gzma",
          "Il1rl1",  #aka St2
          "Il15ra",
          "Rora",
          "S100a4",
          "Tbx21",
          "Tnfrsf25"
          ) 
        ,x.lab.rot = TRUE, point.size.use = 0.5)
dev.off()

VlnPlot(seu,
        c(
          "Gata3",
          "Rora",
          "S100a4",
          "Tbx21"
        ) 
        ,x.lab.rot = TRUE, point.size.use = 0.5)



idtreg <- "Treg"


selcell <- colnames(seu@raw.data)[seu@meta.data$Cd4_exprs_pos & seu@meta.data$doublet<1 & 
                                    seu@meta.data$days.after.infection>=30 & clusterident_mod==idtreg]
#VlnPlot(seu, c("Il10"),x.lab.rot = TRUE, cells.use=selcell)







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

# subsets_cluster <- c("Th1_day_30","Th2_day_30","Treg_day_7_30","Naive_day_7_30")
# subsets_name <- c("Th1","Th2","Treg","Naive")

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




