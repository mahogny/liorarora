#######################################################################
################## Install packages ###################################
#######################################################################
if(FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  biocLite("topGO")
  biocLite("GO.db")
}


#set.seed(1234567890)
library(Seurat)
library(Matrix)
library(biomaRt)
library(stringr)
library(dplyr)
library(limma)
#library(xlsx)

#######################################################################
################## Common GO functions ################################
#######################################################################

library(topGO)
library(GO.db)
require(org.Mm.eg.db)
mapGoTerm <- toTable(GOTERM) 

#####################################
# Human GO database
annotation.mm <- list(
  mapping="org.Mm.eg.db",
  db=org.Mm.eg,
  dbSYMBOL=org.Mm.egSYMBOL,
  dbENSEMBL=org.Mm.egENSEMBL,
  egGO=org.Mm.egGO,
  gene2ez=revmap(org.Mm.egSYMBOL)
)

#####################################
# Get genes for a GO category
go2genes <- function(go,with.descendants=FALSE,annotation){
  genes=c()
  if(with.descendants){
    #Check if recursive
    goids<-go.descendants(go)
  }else{
    goids<-c(go)
  }
  ezgenes<-unique(unlist(mget(goids,revmap(annotation$egGO),ifnotfound=NA)))
  ezgenes<-ezgenes[!is.na(ezgenes)]
  symbols<-unique(unlist(mget(ezgenes,annotation$dbENSEMBL))) #dbSYMBOL
  return(symbols)
}


#####################################
## Function: Get descendant GO categories
go.descendants <- function(go,onto=c("BP","MF","CC")){
  onto=match.arg(onto)
  if(onto=="BP"){
    onto=GOBPOFFSPRING
  }else if(onto=="MF"){
    onto=GOMFOFFSPRING
  }else if(onto=="CC"){
    onto=GOCCOFFSPRING
  }
  children<-unique(unlist(mget(go,onto,ifnotfound=NA)))
  children=union(children,go)
  return(children[complete.cases(children)])
}


stopgo <- function(all_data,DE_data,ID=c("genename","symbol","EnsemblID"),useTopgo=TRUE, keepCol=FALSE,cutoff=1e-4){
  relevant.genes <- rep(1,length(all_data))
  relevant.genes[all_data %in% DE_data] <- 0
  names(relevant.genes) <- all_data
  relevant.genes <- as.factor(relevant.genes)
  
  GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes,nodeSize=5,annot = annFUN.org,
                   mapping="org.Mm.eg.db", ID=ID, geneSel = function(p) p < 0.01)
  # apply Fisher's exact test with elimination mode:
  results <- runTest(GOdata.BP, algorithm = 'elim', statistic = 'fisher')
  x<-GenTable(GOdata.BP, results, topNodes = 60)#, pvalCutOff=0.1)   #was 40
  
  #### TODO: rename things
  colnames(x)[which(colnames(x)=="Term")] <- "term.name"
  colnames(x)[which(colnames(x)=="result1")] <- "p.value"
  x
}


stopgosym <- function(genelist,bg=unique(ensconvert$mgi_symbol),nofactor=FALSE, useTopgo=TRUE,cutoff=1e-4){
  stopgo(bg, genelist, "symbol", useTopgo = useTopgo, cutoff = cutoff) #or genename
}



######################################################################
### Common geneid stuff ##############################################
######################################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]

human_ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#Convert to gene symbols, or retain ID if no symbol
ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart=ensembl),stringsAsFactors=FALSE)

human_ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart=human_ensembl),stringsAsFactors=FALSE)
colnames(human_ensconvert)[2] <- "mgi_symbol"  #should use a different name. genesymbol
v<-sqldf("select *,count(mgi_symbol) as c from human_ensconvert group by ensembl_gene_id")
human_ensconvert <- human_ensconvert[human_ensconvert$ensembl_gene_id %in% v$ensembl_gene_id[v$c==1],] #there are some bastards including CCL3L3. this one better inserted by hand
human_ensconvert$mgi_symbol <- normalizesym(human_ensconvert$mgi_symbol)





nametogenesym<-function(x){
  names(x)<-togenesym(names(x))
  x
}

mtogenesym<-function(x){
  colnames(x)<-togenesym(colnames(x))
  rownames(x)<-togenesym(rownames(x))
  x
}






togenesym <- function(geneid, ensconvert_=ensconvert, dowarn=TRUE){
  ensconvert_ <- ensconvert_[ensconvert_$ensembl_gene_id %in% geneid,]
  ensconvert_ <- ensconvert_[order(ensconvert_$mgi_symbol),]
  ensconvert_ <- ensconvert_[order(ensconvert_$ensembl_gene_id),]
  ensconvert_ <- ensconvert_[!duplicated(ensconvert_$ensembl_gene_id, fromLast=TRUE),]
  rownames(ensconvert_) <- ensconvert_$ensembl_gene_id
  
  out <- ensconvert_[geneid,]$mgi_symbol
  out[is.na(out)] <- geneid[is.na(out)]
  out  
}

toensid <- function(geneid, ensconvert_=ensconvert, dowarn=TRUE){
  ensconvert_ <- ensconvert_[ensconvert_$mgi_symbol %in% geneid,]
  ensconvert_ <- ensconvert_[order(ensconvert_$mgi_symbol),]
  ensconvert_ <- ensconvert_[order(ensconvert_$ensembl_gene_id),]
  ensconvert_ <- ensconvert_[!duplicated(ensconvert_$mgi_symbol, fromLast=TRUE),]
  rownames(ensconvert_) <- ensconvert_$mgi_symbol
  
  out <- ensconvert_[geneid,]$ensembl_gene_id
  out[is.na(out)] <- geneid[is.na(out)]
  out  
}


########################################
## Orthology table human<->mouse
# ortho_mouse_human <- read.csv("mouse_human_ortholog.csv",stringsAsFactors = FALSE)
# colnames(ortho_mouse_human) <- c("ens_mouse","ens_human")
# ortho_mouse_human <- ortho_mouse_human[ortho_mouse_human$ens_human!="" & ortho_mouse_human$ens_mouse!="",]
# ## Only 1-1 mappings
# ortho_mouse_human_unique <- ortho_mouse_human[
#   isUnique(ortho_mouse_human$ens_human) & 
#     isUnique(ortho_mouse_human$ens_mouse),]
# sum(duplicated(ortho_mouse_human_unique$ens_mouse)) 
# sum(duplicated(ortho_mouse_human_unique$ens_human)) 




######################################################################
### Load data ########################################################
######################################################################

######## Load counts and metadata after QC
df <- readRDS('nipo_ss2/df1.rds')
# QC already done

seu<-CreateSeuratObject(raw.data = df$cts, min.cells = 2, project = "T cell", 
                        names.field = 2, names.delim = "\\-")
for(feat in colnames(df$phn)){
  seu <- AddMetaData(seu, setNames(df$phn[[feat]], rownames(df$phn)), feat)
}
seu <- NormalizeData(object = seu)
seu <- FindVariableGenes(object = seu, x.low.cutoff = 0,x.high.cutoff = 5, y.cutoff = 0.1)
length(x = seu@var.genes)


######################################################################
### Annotate particular expression levels ############################
######################################################################
seu <- AddMetaData(seu, setNames(seu@raw.data["Rora",], colnames(seu@raw.data)), "Rora_exprs")
seu <- AddMetaData(seu, setNames(seu@raw.data["Rora",]>5, colnames(seu@raw.data)), "Rora_exprs_pos")
#seu <- AddMetaData(seu, setNames(seu@meta.data$exon4>0, colnames(seu@raw.data)), "exon4_pos")
hist(seu@raw.data["Rora",],breaks=100)
#table(seu@meta.data[,c("treatment","genotype","Rora_exprs_pos")])  #,"exon4_pos"

seu <- AddMetaData(seu, setNames(seu@raw.data["Cd4",]>1, colnames(seu@raw.data)), "Cd4_exprs_pos")
seu <- AddMetaData(seu, setNames(seu@raw.data["Sell",]>1, colnames(seu@raw.data)), "Sell_exprs_pos")
#hist(seu@raw.data["Cd4",],breaks=100)


seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="lungs",  colnames(seu@raw.data)), "is_lung")
seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="spleen", colnames(seu@raw.data)), "is_spleen")
seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="MedLN",  colnames(seu@raw.data)), "is_MedLN")
seu <- AddMetaData(seu, setNames(seu@meta.data$tissue=="MLN",    colnames(seu@raw.data)), "is_MesLN")

# 
# ######################################################################
# ### Only consider CD4+ ###############################################
# ######################################################################
# #seu2 <- SubsetData(object = seu, cells.use = seu@cell.names[seu@meta.data$Cd4_exprs_pos])
# seu2 <- SubsetData(object = seu)#, cells.use = seu@meta.data$Cd4_exprs_pos)
# seu@meta.data$Cd4_exprs_pos
# seu
# seu2 <- SubsetData(object = seu, subset.name = "Cd4_exprs_pos")
# seu2
# 
# seu
# #seu@meta.data$treatment == "uninfected"

######################################################################
### Cluster cells ####################################################
######################################################################
set.seed(1234567890)
seu <- ScaleData(object = seu, genes.use = seu@var.genes, model.use = "negbinom")
seu <- RunPCA(object = seu, pc.genes = seu@var.genes, pcs.compute = 50, pcs.print = 1:2, maxit = 500, weight.by.var = FALSE)
seu <- RunTSNE(object = seu, dims.use = 1:25, do.fast = TRUE)
seu <- FindClusters(object = seu, reduction.type = "pca",  dims.use = 1:25, resolution = 0.5, save.SNN = TRUE, force.recalc = T )
TSNEPlot(object = seu, do.label = TRUE)


#mapcol <- 
#hist(as.double(seu@scale.data[togenesym("Tbx21"),]))

pdf("newclust_allcell.pdf")
TSNEPlot(object = seu, do.label = TRUE)#, colors.use = rep("yellow",1670))
dev.off()

# pdf("newclust_allcell_byday.pdf")
# thecol <- rep("blue",1670)
# thecol[seu@meta.data$days.after.infection==7] <- "red"
# thecol[seu@meta.data$days.after.infection==30] <- "green"
# TSNEPlot(object = seu, colors.use=thecol)
# dev.off()

table(seu@meta.data$days.after.infection)
table(seu@meta.data$tissue)

pdf("newclust_allcell_markers1.pdf",h=20,w=20)
FeaturePlot(seu, c(
  "Tbx21", "Gata3","Il4","Il13","Ifng", #Th1-2 axis
  "Foxp3",       #Tregs
  "Sell","Cd44", #Activated T cells
  "Rorc","Rora",
  #"Il1b",
  #"Ear2",
  
  "Cd19",
  #Cd19, Cd20, Cd34, Cd38, and Cd45r  #immature B cell
  #"Itgax",
  #"Col14a1","Trim47","Ccl12","Ccl7","Apoe","Mgl2","Cd63","Ly86","Cxcl16","Ccl8","Ccl2","Anxa3","Ccl3",
  "Cd68",                         #pulls out macrophages
  #"C1qa","C1qb","C1qc","Sirpa",  #pulls out macrophages
  #"Il17a","Il17f","Il21","Il22",
  #"Ptprc", #Cd45
  #"Cd8a",
  #"Cd14","Lyz2","Ms4a7",
  
  #http://www.abcam.com/primary-antibodies/effector-t-cell-markers
  
  # 	CD8+, CD45RA+, CD45RO-, CCR7+, CD28+	IFN-γ+, IL-2+    Naive Cd8
  #   CD8+, CCR7-	IFN-γ+, Perforin+, Granzyme+    cytotoxic cd8
  # 	CD4+, CCR3+, CCR6+	IL-9+
  #	CD4+, CCR6+, CCR4+, NK1.1+	IL-17+
  #	CD4+, CCR10+, CCR4+, CCR6+ Il22   Th22
  "Il23r",#gammadelta
  "Il21","Bcl6","Cxcr5","Cd40lg","Icos",  #Tfh markers

  "Il1rl1",   #nuocyte markers, these are mixed with the Th2   IL-17BR, ICOS or ST2=il1rl1
  #IL-17A, IL-17F, IL-21, and IL-22
  "days.after.infection", "is_lung","is_spleen","is_MedLN","is_MesLN"), 
  cols.use = c("green", "blue"),pt.size = 1)
dev.off()

#Note: NKT cells can be detected using tracer

# Effector   B cell?

grep("Cd45",ensconvert$mgi_symbol)

## Store cluster identity for later
clusterident <- seu@ident

m <- FindMarkers(seu, ident.1 = "8")
m

#seu

VlnPlot(seu, c("Gata3","Sell","Tbx21","Cd4","Foxp3"))

#TSNEPlot
#?TSNEPlot

#This should contain Th1, Th2, Treg, ... Th17?? 

#seu@scale.data[1:5,1:5]

#8 clusters



######################################################################
### Get markers for each cluster and annotate cell types #############
######################################################################

#seu@raw.data["Rora",]

markers <- FindAllMarkers(object = seu, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers,'outnew_clustermarkers.csv')
#markers = read.csv('markers1.csv')

markers
top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top50 <- markers %>% group_by(cluster) %>% top_n(50, avg_diff)


write.csv(as.data.frame(top20),'outnew_clustermarkers_top20.csv')

#unique(markers$cluster)
#print(head(markers, 5))



######################################################################
### Compare: infected KO CD4+ vs infected control CD4+ ###############
######################################################################


#> CD4+ clusters, treatment: infected, exp: KO, days.after.infection:30, Rora>1, exon4 in KO <1, tissue: lung

### Split cells into 2 groups (and neither)
newident <- rep(0,length(seu@ident))
newident[seu@meta.data$Cd4_exprs_pos & 
         seu@meta.data$days.after.infection>=30 &
         seu@meta.data$genotype=="rora (CD4) KO"] <- "ko"
newident[seu@meta.data$Cd4_exprs_pos & 
         seu@meta.data$days.after.infection>=30 &
         seu@meta.data$genotype=="control mouse"] <- "ctrl"
newident <- as.factor(newident)
sum(newident=="ko")
sum(newident=="ctrl")
names(newident) <- names(seu@ident)
seu@ident <- as.factor(newident)

#as.factor(newident)

### Perform the DE
de_rora <- FindMarkers(seu, ident.1 = "ko", ident.2 = "ctrl")  #, test.use = "DESeq2"

sum(de_rora$p_val<0.05)
de_rora
#Xbp1 is affected with 0.02

write.csv(de_rora, file="out.niposs2/rora_de_genelist.csv")




######################################################################
########### Perform GO analysis
######################################################################
rora_de_go <- stopgosym(
  rownames(de_rora)[de_rora$p_val<0.01],
  rownames(de_rora)
)
#rora_de_go
write.csv(rora_de_go, file="out.niposs2/rora_de_go_p0.01.csv")
rora_de_go <- stopgosym(
  rownames(de_rora)[de_rora$p_val<0.1],
  rownames(de_rora)
)
#rora_de_go
write.csv(rora_de_go, file="out.niposs2/rora_de_go_p0.10.csv")


######################################################################
### Compare: infected KO CD4+ vs infected control CD4+, for each cell type (Th1, Th2, Treg)
######################################################################


#> CD4+ clusters, treatment: infected, exp: KO, days.after.infection:30, Rora>1, exon4 in KO <1, tissue: lung



getdeoneth <- function(theid){
  ### Split cells into 2 groups (and neither)
  newident <- rep(0,length(seu@ident))
  newident[clusterident==theid &
             seu@meta.data$Cd4_exprs_pos & 
             seu@meta.data$days.after.infection>=30 &
             seu@meta.data$genotype=="rora (CD4) KO"] <- "ko"
  newident[clusterident==theid &
             seu@meta.data$Cd4_exprs_pos & 
             seu@meta.data$days.after.infection>=30 &
             seu@meta.data$genotype=="control mouse"] <- "ctrl"
  newident <- as.factor(newident)
  sum(newident=="ko")
  sum(newident=="ctrl")
  names(newident) <- names(seu@ident)
  seu@ident <- as.factor(newident)
  
  ### Perform the DE. for this number ideally with DEseq
  #guess I need to make an object with rounded read counts?
  one_de <- FindMarkers(seu, ident.1 = "ko", ident.2 = "ctrl")  #, test.use = "DESeq2"
  one_de  
}

idth1 <- "4"
idth2 <- "0"
idtreg <- "2"
idnaive <- "1"

de_rora_th1  <- getdeoneth(idth1)
de_rora_th2  <- getdeoneth(idth2)
de_rora_treg <- getdeoneth(idtreg)
de_rora_naive <- getdeoneth(idnaive)   #none - simple!

anydename <- unique(c(rownames(de_rora_th1),rownames(de_rora_th2),rownames(de_rora_treg)))

de_rora_any <- data.frame(
  row.names = anydename, 
  th1p =anydename %in% rownames(de_rora_th1),
  th2p =anydename %in% rownames(de_rora_th2),
  tregp=anydename %in% rownames(de_rora_treg)
)

de_rora_venn <- data.frame(
  row.names = anydename, 
  th1 =anydename %in% rownames(de_rora_th1)[de_rora_th1$p_val<0.01],
  th2 =anydename %in% rownames(de_rora_th2)[de_rora_th2$p_val<0.01],
  treg=anydename %in% rownames(de_rora_treg)[de_rora_treg$p_val<0.01]
)

vc <- vennCounts(de_rora_venn)
vennDiagram(vc,cex=c(1.5,1.5,1.5))

de_rora_venn[apply(de_rora_venn,1,sum)>=2,]
#Etf1, H2afz, Id2

##################### compare fractions, rora/ctrl #######################

isrel <-
 clusterident %in% c(idth1,idth2,idtreg) &  ########## todo: add naive to this list!
  seu@meta.data$Sell_exprs_pos &
  seu@meta.data$Cd4_exprs_pos & 
  seu@meta.data$Cd4_exprs_pos & 
  seu@meta.data$days.after.infection>=30
#seu@meta.data$genotype %in% c("control mouse", "rora (CD4) KO")

mousenum_ctrl <- unique(seu@meta.data$mouse.number[seu@meta.data$genotype=="control mouse"])
mousenum_ko   <- unique(seu@meta.data$mouse.number[seu@meta.data$genotype=="rora (CD4) KO"])
 
compnum <- table(seu@meta.data$mouse.number[isrel], clusterident[isrel])
compnum <- compnum[,c(idth1,idth2,idtreg)] 
for(i in 1:nrow(compnum)){
  compnum[i, ] <- compnum[i, ]/sum(compnum[i, ])
}
compnum_ctrl <- compnum[rownames(compnum) %in% mousenum_ctrl,]
compnum_ko   <- compnum[rownames(compnum) %in% mousenum_ko,]

compnum_ctrl_mean <- apply(compnum_ctrl,2,mean)
compnum_ko_mean   <- apply(compnum_ko,  2,mean)   #KO has a slight increase in Th2/Treg than control

compnum_ctrl_mean[2]/compnum_ctrl_mean[3]
compnum_ko_mean[2]/compnum_ko_mean[3]

t.test(compnum_ctrl[,2]/compnum_ctrl[,3], compnum_ko[,2]/compnum_ko[,3])
#t.test(compnum_ctrl[,3], compnum_ko[,3])


########### what about the number of activated cells????




isrel <-
  clusterident %in% c(idth1,idth2,idtreg) &
#  seu@meta.data$Sell_exprs_pos &
  seu@meta.data$Cd4_exprs_pos & 
  seu@meta.data$Cd4_exprs_pos & 
  seu@meta.data$days.after.infection>=30
#seu@meta.data$genotype %in% c("control mouse", "rora (CD4) KO")

mousenum_ctrl <- unique(seu@meta.data$mouse.number[seu@meta.data$genotype=="control mouse"])
mousenum_ko   <- unique(seu@meta.data$mouse.number[seu@meta.data$genotype=="rora (CD4) KO"])

compnum <- table(seu@meta.data$mouse.number[isrel], clusterident[isrel])
compnum <- compnum[,c(idth1,idth2,idtreg)] 
for(i in 1:nrow(compnum)){
  compnum[i, ] <- compnum[i, ]/sum(compnum[i, ])
}
compnum_ctrl <- compnum[rownames(compnum) %in% mousenum_ctrl,]
compnum_ko   <- compnum[rownames(compnum) %in% mousenum_ko,]

compnum_ctrl_mean <- apply(compnum_ctrl,2,mean)
compnum_ko_mean   <- apply(compnum_ko,  2,mean)   #KO has a slight increase in Th2/Treg than control

compnum_ctrl_mean[2]/compnum_ctrl_mean[3]
compnum_ko_mean[2]/compnum_ko_mean[3]



compnum[mousenum_ctrl,]

compnum

mousenum_ko

 #?table
 
unique(seu@meta.data$mouse.number)

seu@meta.data$

#how many 



######################################################################
### Compare with RORA motifs, any particular direct targets? 
######################################################################

motifu <- read.csv("out/upstream.csv",sep="\t",stringsAsFactors = FALSE)  
motifd <- read.csv("out/downstream.csv",sep="\t",stringsAsFactors = FALSE)  

#5kb analyzed. subset here
keepdist <- 5000
motifud <- rbind(
  motifu[motifu$start>5000-keepdist,],
  motifu[motifu$end<keepdist,])

motifud$gene <- str_split_fixed(motifud$sequence_name,":",2)[,1]
motifud$gene

#Extract one motif
#usemotif <- "MA0071.1"
usemotif <- "MA0072.1"
genesWithMotif <- motifud$gene[motifud$X..motif_id=="MA0072.1" & motifud$p.value<1e-5]
genesWithMotifSym <- ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% genesWithMotif]

length(genesWithMotif)  #5692
sum(de_rora$p_val<0.01) #41

intersect(rownames(de_rora)[de_rora$p_val<0.01], genesWithMotifSym)

#S100a10"  "S100a4"   "Aff4"     "Slc25a19" "Npm1"
#in vitro: aff4 has the same trend as rora, not really different Th0/2
#in vitro: s100a4 decreases, different Th0/2
#in vitro: Slc25a19 similar to rora, but does not go to 0
#in vitro: s100a10 increases continuously, different Th0/2
#in vitro: Npm1 goes up toward d3, different Th0/2
#s100a4 is a strong hit for gata3


