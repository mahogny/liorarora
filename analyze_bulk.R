#######################################################################
################## Install packages ###################################
#######################################################################

library(stringr)
library(limma)
library(sqldf)
library(DESeq2)
library(DESeq)

source("common_functions.R")

ens_il10 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Il10"]
ens_mboat1 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Mboat1"]


#######################################################################
################## Analyze RORA siRNA in human Th17 ###################
#######################################################################

read_th17 <- function(){
  dat <- read.csv("ext_data/th17_count.csv", sep="\t", stringsAsFactors = FALSE)
  cond <- read.csv("ext_data/th17_cond.csv", sep="\t", stringsAsFactors = FALSE)
  cond_libname <- read.csv("ext_data/th17_cond_libraryname.csv", stringsAsFactors = FALSE)
  cond_libname$Sample_Name <- str_trim(cond_libname$Sample_Name)
  cond <- merge(cond,cond_libname)
  rownames(cond) <- cond$libraryname
  cond <- cond[order(rownames(cond)),]
  dat <- dat[,order(colnames(dat))]
  rownames(dat) <- dat$geneid
  dat <- dat[,colnames(dat) %in% rownames(cond)]
  list(count=dat, cond=cond)
}

dat_th17 <- read_th17()

#Locate relevant conditions
dat_th17$cond$rora_treatment <- NA
dat_th17$cond$rora_treatment[dat_th17$cond$source_name=="CD4 T cell_Mock_24h"] <- "ctrl"
dat_th17$cond$rora_treatment[dat_th17$cond$treatment=="RORA siRNA_1"] <- "ko"
dat_th17$cond$rora_treatment[dat_th17$cond$treatment=="RORA siRNA_2"] <- "ko"
dat_th17$cond$rora_treatment <- factor(dat_th17$cond$rora_treatment)

#Find DE genes for the siRNA
keep <- !is.na(dat_th17$cond$rora_treatment)
dds <- DESeqDataSetFromMatrix(countData = round(dat_th17$count[,keep]),
                              colData = dat_th17$cond[keep,],
                              design = ~  rora_treatment)
dds <- DESeq(dds)
v <- results(dds)
v <- as.data.frame(v[order(v$pvalue),])

#Extract human gene name
v$mgi_symbol <- normalizesym(str_split_fixed(rownames(v),"_",2)[,1])

#Convert to gene ID
v <- merge(v, human_ensconvert)
colnames(v)[colnames(v)=="ensembl_gene_id"] <- "ens_human"

#Convert to mouse gene ID
v <- merge(v, ortho_mouse_human)
#colnames(v)[colnames(v)=="ens_mouse"] <- "ensembl_gene_id"
#v <- merge(v,ensconvert)
#v$mgi_symbol

#v$ens_human
de_th17 <- data.frame(
  p_th17=v$pvalue, 
  fc_th17=v$log2FoldChange, 
  ensembl_gene_id=v$ens_mouse,
  stringsAsFactors = FALSE)



########## how similar are the DE gene list, rora1 vs rora2?



#######################################################################
############ Analyze Rora KO in liver macrophages #####################
#######################################################################

read_liverko <- function(){
  dat <- read.csv("ext_data/fpkm.liver_roraKO.csv", sep="\t", stringsAsFactors = FALSE)
  colnames(dat)
  rownames(dat) <- dat$gene_id
  dat <- dat[,-1]
  
  cond <- data.frame(treatment=rep("ctrl",8), stringsAsFactors = FALSE)
  cond$treatment[grep("KO",colnames(dat))] <- "ko"
  cond$treatment <- factor(cond$treatment)
  list(count=dat, cond=cond)
}
dat_liverko <- read_liverko()

#Find DE genes
dds <- DESeqDataSetFromMatrix(countData = round(dat_liverko$count),
                              colData = dat_liverko$cond,
                              design = ~  treatment)
dds <- DESeq(dds)
v <- as.data.frame(results(dds))
v <- v[order(v$pvalue),]
de_livermacro <- data.frame(
  p_macro=v$pvalue, 
  fc_macro=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),
  stringsAsFactors = FALSE)


#head(de_livermacro)


######################################################################
### Overexpression data ##############################################
######################################################################

#### Read count matrix
datOver <- read.csv("overexpression/roraOverexprCount.csv",sep="\t")
rownames(datOver) <- datOver[,1]
datOver <- datOver[,-1]
datOver <- datOver[-grep("_",rownames(datOver)),]  ##remove __no_feature, etc

#### Build cell condition matrix
datOverCond <- str_split_fixed(colnames(datOver),"_",2)
colnames(datOverCond) <- c("v","condition")
datOverCond <- as.data.frame(datOverCond)
datOverCond$mouse <- str_sub(datOverCond$v,2,2)
datOverCond$mouse <- as.factor(datOverCond$mouse)

write.csv(datOverCond, "paper//supplementary files/s2_rnaseq/bulk_overexpression/cellcond.csv")

#M6 GFP is the reference level for condition
dds <- DESeqDataSetFromMatrix(countData = datOver,
                              colData = datOverCond,
                              design = ~condition + mouse)
dds <- DESeq(dds)

getroraoverexpDE <- function(contr){
  ### FC is M6_GFP / contr
  v <- results(dds, contrast = c("condition","M6_GFP",contr)) ########## wait. why not use M6_mCh as well?
  v <- v[order(v$pvalue),]
  v$pvalue[is.na(v$pvalue)] <- 1
  v  
}
roratopDE <- function(v){
  v <- as.data.frame(v[1:100,])
  v$genesym <- togenesym(rownames(v))
  v  
}


#### Calculate results for each construct. V1 and V3 agrees the most according to venn diagram
v1 <- getroraoverexpDE("rora1_GFP")
v2 <- getroraoverexpDE("rora1_mCh")  #strong p-values. Cd84, Cd27, Ccr5, Cd63, Il4 etc  ... see which one overlaps best?
v3 <- getroraoverexpDE("rora2_mCh")  #rora is changed!  il10. 

#### Check how the different DE compare
pcut <- 1e-5
allu <- unique(c(rownames(v1), rownames(v2), rownames(v3)))
# overexpvenn <- data.frame(
#   row.names = allu, 
#   v1=allu %in% rownames(v1)[v1$pvalue<pcut],
#   v2=allu %in% rownames(v2)[v2$pvalue<pcut],
#   v3=allu %in% rownames(v3)[v3$pvalue<pcut]
# )
overexpvenn <- data.frame(
  row.names = allu, 
  v1=allu %in% rownames(v1)[rank(v1$pvalue)<100],   ###Base it on FC instead
  v2=allu %in% rownames(v2)[rank(v2$pvalue)<100],
  v3=allu %in% rownames(v3)[rank(v3$pvalue)<100]
)
vc <- vennCounts(overexpvenn)
vennDiagram(vc,cex=c(1.5,1.5,1.5))

#### Genes to consider later for overlapping
overexp_result <- roratopDE(v3)
v <- v3
de_overexp <- data.frame(
  p_overexp=v$pvalue, 
  fc_overexp=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),
  stringsAsFactors = FALSE)
#roratopDE(v3)

#v3 is n=2+3
pdf(sprintf("out.de_comp/newvolc_oe.pdf"),width = 4,height = 5)
new_volcano(v[v$baseMean>20,])
dev.off()
#

v[rownames(v)==toensid("Id2"),]
v[rownames(v)==toensid("Myc"),]
v[rownames(v)==toensid("St6galnac6"),]
v[rownames(v)==toensid("Rora"),]  #Rora up => Id2 and Myc up, and st6..

togenesym(rownames(v[v$baseMean>50,])[1:20])   
#Serinc2, CCdc184, St6galnac6 !!!!, Ccr5, Myc, Havcr2



########## test limma here!!!!


# 
# library(edgeR)
# #keepdat <- datOverCond$condition %in% c("M6_GFP","rora1_GFP")
# keepdat <- datOverCond$condition %in% c("M6_GFP","M6_mCh","rora1_mCh","rora2_mCh")
# newcond <- datOverCond[keepdat,]
# newcond$isko <- 0
# newcond$isko[newcond$condition %in% c("rora1_mCh","rora2_mCh")] <- 1
# newcond
# newcond$condition <- factor(as.character(newcond$condition))
# dge <- DGEList(counts=datOver[,keepdat])
# des <- model.matrix(data=newcond, ~condition )  #+ mouse
# keep <- filterByExpr(dge, des)
# dge <- dge[keep,,keep.lib.sizes=FALSE]
# dge <- calcNormFactors(dge)
# logCPM <- cpm(dge, log=TRUE, prior.count=3)
# fit <- lmFit(logCPM, des)
# fit <- eBayes(fit, trend=TRUE)
# v <- topTable(fit, coef=ncol(des), number = 1000)
# v$genesym <- togenesym(rownames(v))
# v

#not in any with ~condition + mouse.
#is in rora2_mCh, no mouse
"ENSMUSG00000024793" %in% rownames(v)    #Tnfrsf25
#also in top 1000, rora1_mCh

                       
                       
#########################################################################
################## Assemble condition matrix, 41 ########################
#########################################################################


cellcond41 <- read.csv("walkup/walkup41/cellcond41.csv",stringsAsFactors = FALSE, row.names = 1)
cellcond41[is.na(cellcond41)] <- ""
cellcond41 <- cellcond41[order(cellcond41$index),]  

cellcond41$altko <- factor(cellcond41$altko)
cellcond41$isko <- factor(cellcond41$isko)
cellcond41$isinf <- factor(cellcond41$isinf)
cellcond41$mouse <- factor(cellcond41$mouse)
cellcond41$isnewspleen <- factor(cellcond41$isnewspleen)
#alltissues <- unique(cellcond41$tissue[cellcond41$tissue!=""])


#########################################################################
###################### Read counts and do QC, 41 ########################
#########################################################################

## Read exon4
#cnt41ex <- read.csv("walkup/counttable_exon.csv", row.names = 1)
#cellcond41$exon4 <- cnt41ex[grep("ENSMUSE00001263945",rownames(cnt41ex)),]

## Read counts
if(FALSE){
  cnt41 <- read.csv("walkup/walkup41/counttable.csv", row.names = 1)
  cnt41 <- cnt41[,colnames(cnt41) %in% rownames(cellcond41)]
  #rownames(cellcond41) <- colnames(cnt41)   
  nrow(cellcond41)
  
  ## Put transcripts together into genes
  cnt41 <- rowsum(cnt41, ens_transc_gene[str_split_fixed(rownames(cnt41), pattern = "\\.", 2)[,1],]$ensembl_gene_id)
  
  apply(cnt41,2,sum)  ###highly suspicious
  
  #Filter out lowly expressed genes
  sum(apply(cnt41,1,mean)>5)
  cnt41 <- cnt41[apply(cnt41,1,mean)>5,]
  #sum(apply(cnt41>5,1,any))
  #sum(apply(cnt41>0,1,any))
} else {
  
  cnt41 <- read.csv("walkup/walkup41/counttable_ht.csv", row.names = 1)
  cnt41 <- cnt41[,colnames(cnt41) %in% rownames(cellcond41)]
  
  #Filter out lowly expressed genes
  sum(apply(cnt41,1,mean)>5)
  cnt41 <- cnt41[apply(cnt41,1,mean)>5,]
}


write.csv(cellcond41, "out_cellcond41.csv")


#########################################################################   all activated cells?
############################### Nipo day 30 #############################   at what date?
#########################################################################


keep <- which(cellcond41$isnewspleen=="0" & cellcond41$isinf=="1")
dds <- DESeqDataSetFromMatrix(countData = round(cnt41[,keep]),
                              colData = cellcond41[keep,c("isko"),drop=FALSE],
                              design = ~ isko)
dds <- DESeq(dds)
v <- results(dds)
v <- v[order(v$pvalue),]

de_nipo30 <- data.frame(
  p_nipo30=v$pvalue, 
  fc_nipo30=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),stringsAsFactors = FALSE)



#########################################################################   all activated cells?
############################### Nipo day 0 ##############################    at what date?
#########################################################################


keep <- which(cellcond41$isnewspleen=="0" & cellcond41$isinf=="0")
dds <- DESeqDataSetFromMatrix(countData = round(cnt41[,keep]),
                              colData = cellcond41[keep,c("isko"),drop=FALSE],
                              design = ~ isko)
dds <- DESeq(dds)
v <- results(dds)
v <- v[order(v$pvalue),]

de_nipo0 <- data.frame(
  p_nipo0=v$pvalue, 
  fc_nipo0=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),stringsAsFactors = FALSE)



#########################################################################  
############################### Tregs only, day 0, spleen ###############   
#########################################################################


keep <- which(cellcond41$isnewspleen=="1" & cellcond41$isko %in% c("0","1"))
dds <- DESeqDataSetFromMatrix(countData = round(cnt41[,keep]),
                              colData = cellcond41[keep,c("isko"),drop=FALSE],
                              design = ~ isko)
dds <- DESeq(dds)
v <- results(dds)
v <- v[order(v$pvalue),]

de_treg_spl_0 <- data.frame(
  p_treg_spl_0=v$pvalue, 
  fc_treg_spl_0=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),stringsAsFactors = FALSE)





#########################################################################
################## Assemble condition matrix, 37 ########################
#########################################################################

cellcond37 <- read.csv("walkup/walkup37/cellcond.csv",stringsAsFactors = FALSE, row.names = 1)
cellcond37$isko <- factor(cellcond37$isko)
cellcond37$group <- factor(cellcond37$group)
cellcond37$tissue <- factor(cellcond37$tissue)
cellcond37$mouse <- factor(cellcond37$mouse)
alltissues <- unique(cellcond37$tissue[cellcond37$tissue!=""])


#########################################################################
###################### Read counts and do QC, 37 ########################
#########################################################################

if(FALSE){
  
  ## Read counts
  cnt37 <- read.csv("walkup/counttable2.csv", row.names = 1)
  cnt37 <- cnt37[,colnames(cnt37) %in% rownames(cellcond37)]
  #rownames(cellcond37) <- colnames(cnt37)   
  #nrow(cellcond37)
  
  
  ## No idea why NA?
  # keep <- !is.na(cnt37[10,])
  # cellcond37$name[!keep]
  # cnt37 <- cnt37[,keep] 
  # cellcond37 <- cellcond37[keep,]
  
  
  #cnt37 <- round(cnt37)
  
  ## Only with counts > 0
  keep <- apply(cnt37,2,sum)>0    #where did the NA come from?
  cellcond37$name[!keep]  #""          ""          ""          "WT_1"      "WT_4"      "D1 S3 114" "D1 S1 116" ""  
  cnt37 <- cnt37[,keep] 
  cellcond37 <- cellcond37[keep,]
  
  
  ## Put transcripts together into genes
  cnt37 <- rowsum(cnt37, ens_transc_gene[str_split_fixed(rownames(cnt37), pattern = "\\.", 2)[,1],]$ensembl_gene_id)
  
  #Filter out lowly expressed genes
  cnt37 <- cnt37[apply(cnt37,1,mean)>5,]
  nrow(cnt37)
  
} else {
  
  cnt37 <- read.csv("walkup/walkup37/counttable_ht.csv", row.names = 1)
  cnt37 <- cnt37[,colnames(cnt37) %in% rownames(cellcond37)]
  
  ## Only with counts > 0
  keep <- apply(cnt37,2,sum)>0    #where did the NA come from?
  cellcond37$name[!keep]  #"WT_1"      "WT_4"      "D1 S3 114" "D1 S1 116" ""  
  cnt37 <- cnt37[,keep] 
  cellcond37 <- cellcond37[keep,]
  
  #Filter out lowly expressed genes
  sum(apply(cnt37,1,mean)>5)
  cnt37 <- cnt37[apply(cnt37,1,mean)>5,]
}

write.csv(cellcond37, "out_cellcond37.csv")

#########################################################################
##################### Nipo day 7. activated t cells #####################
#########################################################################

keep <- which(cellcond37$group=="oldko")
cellcond37[keep,]
dds <- DESeqDataSetFromMatrix(countData = round(cnt37[,keep]),
                              colData = cellcond37[keep,],
                              design = ~ isko)
dds <- DESeq(dds)
v <- results(dds)
v <- v[order(v$pvalue),]

de_nipo_7 <- data.frame(
  p_nipo_7=v$pvalue, 
  fc_nipo_7=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),stringsAsFactors = FALSE)



#########################################################################
########### Tissues, one by one, day 30, activated CD4+CD62L-CD44+ ######
#########################################################################

allp <- list()
for(i in 1:length(alltissues)){
  keep <- which(cellcond37$group=="newko" & as.character(cellcond37$tissue)==as.character(alltissues[i]))
  dds <- DESeqDataSetFromMatrix(countData = round(cnt37[,keep]),
                                colData = cellcond37[keep,],
                                design = ~ isko )
  dds <- DESeq(dds)
  v <- results(dds, c("isko","1","0"))
  v1 <- data.frame(ensembl_gene_id=rownames(v),stringsAsFactors = FALSE)
  v1[sprintf("p_%s",alltissues[i])] <- v$pvalue
  v1[sprintf("fc_%s",alltissues[i])] <- v$log2FoldChange
  if(i==1){
    allp <- v1
  } else {
    allp <- merge(allp, v1, all=TRUE)
  }
}
allp[is.na(allp)]<-1
rownames(allp) <- allp$trans

de_tissue <- allp ########## check

##############
new_volcano <- function(v){
  to_highlight <- c("Il17a","Il10","Il5","Ccl3","Il12rb2","Cxcr6","St6galnac6","Tnfrsf25","Arntl","S100a4","Alox8","Sell",
                    "Myc","Id2","Cxcr6","Tbx21","Pparg","Fasl","Gzma",
                    "Colq","Ccr6","Emb","Gzmb")  #Emb is relevant!
  
  highlight_g <- c(to_highlight,na.omit(rownames(v)[v$pvalue<1e-20]))
  tohi_ens <- ensconvert[ensconvert$mgi_symbol %in% highlight_g,]$ensembl_gene_id
  
  #Il17a, Il5, Il2 down, while Ccl5, Il12rb2, Ccl3 are up. Furthermore Tbx21
  
  #On top: Hlf. then gzma
  
  v<-as.data.frame(v)
  rem <- c(
    ensconvert[grep("Igh",ensconvert$mgi_symbol),]$ensembl_gene_id,
    ensconvert[grep("Igk",ensconvert$mgi_symbol),]$ensembl_gene_id)
  v<-v[!(rownames(v) %in% rem),]
  v<-v[order(v$pvalue),]
  
  xlim <- range(v$log2FoldChange,na.rm = TRUE)
  ylim <- range(-log10(v$pvalue),na.rm = TRUE)
  
  v3 <- v[!(rownames(v) %in% tohi_ens),]
  plot(v3$log2FoldChange,-log10(v3$pvalue),pch=19,cex=0.3,col="darkgray",
       ylab="Log10 p-value",xlab="Log2 Fold change",
       xlim=xlim, ylim=ylim)
  
  v2 <- v[rownames(v) %in% tohi_ens,]    
  text(x=v2$log2FoldChange, y=-log10(v2$pvalue), labels = togenesym(rownames(v2)))
  
#  v2 <- v[rownames(v) %in% highlight_g,]    
#  text(x=v2$log2FoldChange, y=-log10(v2$pvalue), labels = rownames(v2))
}


if(FALSE){
  #Colon is weird. cytokines coming up. this is a perturbation analysis

  for(i in 1:4){
 #   i <- 3
    print(alltissues[i])
    keep <- which(
        cellcond37$group=="newko" & 
        as.character(cellcond37$tissue)==as.character(alltissues[i]) &
        rownames(cellcond37)!="samp_219" &  #from lung
        rownames(cellcond37)!="samp_267__")
    dds <- DESeqDataSetFromMatrix(countData = round(cnt37[-grep("_",rownames(cnt37)),keep]),  #-grep("_",rownames(cnt37))
                                  colData = cellcond37[keep,],
                                  design = ~ isko )
    dds <- DESeq(dds)
    v <- results(dds, c("isko","1","0"))
    
    ### for i=4, colon, check samp_217
    # counts(dds, normalized=TRUE)[toensid("Tbx21"),]    #super high first sample   217   258
    # counts(dds, normalized=TRUE)[toensid("Cxcr6"),]    #super high first sample
    # counts(dds, normalized=TRUE)[toensid("Il5"),]      #fairly spread
    # counts(dds, normalized=TRUE)[toensid("Tnfrsf25"),]      #fairly spread
    # counts(dds, normalized=TRUE)[toensid("Il10"),]
    # counts(dds, normalized=TRUE)[toensid("Il17a"),]   #why reported??
    # counts(dds, normalized=TRUE)[toensid("Il2"),]   #0 first sample
    # counts(dds, normalized=TRUE)[toensid("Il12rb2"),]  #super high first sample

    # counts(dds, normalized=TRUE)[toensid("Ccnd2"),]  #219 wtf
    # counts(dds, normalized=TRUE)[toensid("Lck"),]  #219 wtf
    # counts(dds, normalized=TRUE)[toensid("Pparg"),]   #219 wtf
    # counts(dds, normalized=TRUE)[toensid("Gnai3"),]   #  267?
    # counts(dds, normalized=TRUE)[toensid("Narf"),]   #  267?
    # counts(dds, normalized=TRUE)[toensid("Itgb2"),]   #  267?
    # counts(dds, normalized=TRUE)[toensid("Sept1"),]   #  
    
    # Il17a, Il5, Il2 down, while Ccl5, Il12rb2, Ccl3 
    # v[toensid("Tbx21"),]  
    # v[toensid("Il10"),]  
    # v[toensid("Il5"),]  
    # v[toensid("Il17a"),]  
    # v[toensid("Ccnd2"),]  
    # tail(cnt37)
    
    print(nrow(cellcond37[keep,]))
    pdf(sprintf("out.de_comp/newvolc_%s.pdf",alltissues[i]),width = 4,height = 5)
    new_volcano(v[v$baseMean>20,])
#    togenesym(rownames(v[v$baseMean>20,])[1:30])
#    v[toensid("Cxcr6"),]
      #Cxcr6! 
    
    dev.off()
        
  }
  
}


# allp <- allp[order(apply(allp[,c("p_spleen","p_gut","p_lung","p_colon")],1,min)),]
# 
# allp$num_sig <- apply(allp[,c("p_spleen","p_gut","p_lung","p_colon")]<0.01,1,sum)
# a <- allp[allp$num_sig>=2,]
# a$genesym <- togenesym(a$ensembl_gene_id)
# a[order(a$num_sig,decreasing = TRUE),]
# #a[,c("genesym","num_sig",colnames(a)[grep("p_",colnames(a))])]
# 
# allp$num_sig <- apply(allp[,c("p_spleen","p_gut","p_lung","p_colon")]<0.05,1,sum)   #0.01
# a <- allp[allp$num_sig>=2,]
# a$genesym <- togenesym(a$ensembl_gene_id)
# #a[order(a$num_sig,decreasing = TRUE),]
# a[order(apply(a[,c("fc_spleen","fc_gut","fc_lung","fc_colon")],1,function(x) mean(na.omit(x))),decreasing = TRUE),]
# 

# 
# 
# ## Venn diagram on gene level
# vc <- vennCounts(allp[,c("p_spleen","p_gut","p_lung","p_colon")]<0.01)
# vennDiagram(vc,cex=c(1.5,1.5,1.5))
# 
# allp$mean_fc <- apply(allp[,grep("fc_",colnames(allp))],1,mean)
# allp <- allp[order(abs(allp$mean_fc),decreasing = TRUE),]
# a <- allp[1:100,]
# a$genesym <- togenesym(a$ensembl_gene_id)
# a
# 
# allp$num_sig <- apply(allp[,c("p_spleen","p_gut","p_lung","p_colon")]<0.01,1,sum)
# a <- allp[allp$num_sig>=2,]
# a$genesym <- togenesym(a$ensembl_gene_id)
# a
# a[,c("genesym","num_sig",colnames(a)[grep("fc_",colnames(a))])]



#########################################################################
########### Tissues, COMBINED, day 30, activated CD4+CD62L-CD44+ ########
#########################################################################

#"spleen" drives DE of gene CD4 ... ,"gut", "lung"
#colon drives Igkv***. and lung. and gut
keep <- which(cellcond37$group=="newko" & cellcond37$tissue %in% c("gut","lung","colon","spleen")) 
cellcond37[keep,]
dds <- DESeqDataSetFromMatrix(countData = round(cnt37[,keep]),
                              colData = cellcond37[keep,],
                              design = ~   isko)
dds <- DESeq(dds)
v <- results(dds, contrast = c("isko","0","1"))
v <- v[order(v$pvalue),]

de_tissue_combined <- data.frame(
  p_tcomb=v$pvalue, 
  fc_tcomb=v$log2FoldChange, 
  ensembl_gene_id=rownames(v),stringsAsFactors = FALSE)





#########################################################################
################## Sci.Imm skin data ####################################
#########################################################################


sci_cellcond <- read.csv("ext_skintreg/skincond.csv")
rownames(sci_cellcond) <- sci_cellcond$name
sci_cellcond$isko <- factor(sci_cellcond$isko)
sci_count <- read.csv("ext_skintreg/GSE99086_Read_Counts_All_Cell_Types.csv", stringsAsFactors = FALSE,sep="\t")
rownames(sci_count) <- sci_count[,1]
sci_count<- sci_count[,-1]
rownames(sci_cellcond) <- colnames(sci_count)


########## Compare Treg wt/KO in tissue, from the other paper
keep <- which(sci_cellcond$cell=="Treg" & sci_cellcond$tissue=="skin")
sci_dds <- DESeqDataSetFromMatrix(countData = sci_count[,keep],
                                  colData = sci_cellcond[keep,],
                                  design = ~ isko) 
sci_dds <- DESeq(sci_dds)
v <- results(sci_dds)

de_skin_treg <- data.frame(
  p_skin_treg=v$pvalue, 
  fc_skin_treg=v$log2FoldChange, 
  mgi_symbol=rownames(v),stringsAsFactors = FALSE)

de_skin_treg_ens <- merge(de_skin_treg, ensconvert)

pdf(sprintf("out.de_comp/newvolc_skin.pdf"),width = 4,height = 5)
new_volcano_skin(v[v$baseMean>20,])
dev.off()


new_volcano_skin <- function(v){
  to_highlight <- c("Il17a","Il10","Il5","Ccl3","Il12rb2","Cxcr6","St6galnac6","Tnfrsf25","Arntl","S100a4","Alox8","Sell",
                    "Myc","Id2","Cxcr6","Tbx21","Pparg","Emb","Gzma","Fasl")
  
  highlight_g <- c(to_highlight,na.omit(rownames(v)[v$pvalue<1e-40]))
  tohi_ens <- ensconvert[ensconvert$mgi_symbol %in% highlight_g,]$ensembl_gene_id
  
  v<-as.data.frame(v)
#  v<-v[rownames(v) %in% ensconvert[-grep("Igh",ensconvert$mgi_symbol),]$ensembl_gene_id,] #contamination spotted
 # v<-v[rownames(v) %in% ensconvert[-grep("Igk",ensconvert$mgi_symbol),]$ensembl_gene_id,] #contamination spotted
  v<-v[order(v$pvalue),]
  
  xlim <- range(v$log2FoldChange,na.rm = TRUE)
  ylim <- range(-log10(v$pvalue),na.rm = TRUE)
  
  v3 <- v[!(rownames(v) %in% highlight_g),]
  plot(v3$log2FoldChange,-log10(v3$pvalue),pch=19,cex=0.3,col="darkgray",
       ylab="Log10 p-value",xlab="Log2 Fold change",
       xlim=xlim, ylim=ylim)
  v2 <- v[rownames(v) %in% highlight_g,]    
  text(x=v2$log2FoldChange, y=-log10(v2$pvalue), labels = rownames(v2))
}


#######################################################################
################## Read RORA KO in human macrophages ##################
#######################################################################



dat_thp1 <- read.csv("crispr_thp1/de_thp1.csv")
dat_thp1$Gene <- normalizesym(dat_thp1$Gene)
colnames(dat_thp1)[1] <- "mgi_symbol"
head(dat_thp1)
head(human_ensconvert)

human_ensconvert$mgi_symbol

dat_thp1 <- merge(dat_thp1, human_ensconvert)
dat_thp1 <- dat_thp1[dat_thp1$mgi_symbol!="",]
dat_thp1 <- sqldf("select distinct * from dat_thp1")
colnames(dat_thp1)[7] <- c("ens_human")

dat_thp1 <- merge(dat_thp1, ortho_mouse_human_unique)


de_thp1 <- data.frame(
  p_thp1=dat_thp1$pval, 
  fc_thp1=dat_thp1$logFC, 
  ensembl_gene_id=dat_thp1$ens_mouse,
  stringsAsFactors = FALSE)

de_thp1
