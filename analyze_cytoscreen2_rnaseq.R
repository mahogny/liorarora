library(gplots)
library(reshape2)
library(scater)
library(stringr)
library(sqldf)
library(DESeq)
library(DESeq2)
library(Rtsne)
library(limma)

source("common_functions.R")

######################################################################
############# Read plate layout ######################################
######################################################################

#Read layout. Correct the layout
cytolay_f <- read.table("upstream2/layout_rnaseq_cytoscreen.csv", 
                        stringsAsFactors = FALSE, sep = ",")
cytolay_mouse <- read.table("upstream2/layout_rnaseq_cytoscreen_mouse.csv", 
                            stringsAsFactors = FALSE, sep = ",")
cytolay <- NULL
for(i in 1:12){
  for(j in 1:8){
    cytolay <- rbind(cytolay,
                     data.frame(well=sprintf("%s%s",c("A","B","C","D","E","F","G","H")[j],i), 
                                treatment=cytolay_f[j,i], 
                                mouse=cytolay_mouse[j,i],
                                stringsAsFactors = FALSE))
  }
}


wellname <- matrix(nrow=8,ncol=12)
for(i in 1:12){
  for(j in 1:8){
    wellname[j,i]<-sprintf("%s%s",c("A","B","C","D","E","F","G","H")[j],i)
  }
}


####### Build up index numbers for rnaseq plates

rnaseqPlateIndex <- read.table("upstream/rnaseqplateindex.csv",stringsAsFactors = FALSE, sep = ",")
rnaseqPlateIndex <- as.matrix(rnaseqPlateIndex)

mapLibindexWell <- data.frame(well=wellname[1:96], libraryid=rnaseqPlateIndex[1:96], plate=1)
mapLibindexWell$libraryname <- sprintf("samp_%s",mapLibindexWell$libraryid)

cellcond <- merge(cytolay, mapLibindexWell)
cellcond <- cellcond[order(cellcond$libraryid),]



######################################################################
############# Read count data and make condition matrix ##############
######################################################################

dat <- read.csv("upstream2/counttable.csv",stringsAsFactors = FALSE, sep = ",", row.names = 1)
rownames(dat) <- str_split_fixed(rownames(dat),"\\.",2)[,1]
dat <- dat[,1:96]

mapGeneidTransid <- read.csv("upstream/map_ensg_ensmust.csv", stringsAsFactors = FALSE)
rownames(mapGeneidTransid) <- mapGeneidTransid$transid
dat_geneid <- mapGeneidTransid[rownames(dat),]$ensemble_gene_id

#Sum on per-gene level. Drop unknown transcripts
dat_geneid[is.na(dat_geneid)] <- "NA"
sum(dat_geneid=="NA")
dat_pergene <- rowsum(dat, group=dat_geneid)
dat_pergene <- dat_pergene[rownames(dat_pergene)!="NA",]
dat <- dat_pergene

dim(dat)
dim(cellcond)

######################################################################
############# Read distribution across plate #########################
######################################################################


##Read distribution across plate
cellcond$total_reads <- apply(dat,2,sum)
outmat <- matrix(nrow=8, ncol=12)
for(i in 1:12){
  for(j in 1:8){
    outmat[j,i] <- cellcond$total_reads[cellcond$well==wellname[j,i]]
  }
}
round(outmat)





######################################################################
############# QC #####################################################
######################################################################


#exclude wells known to be empty
keep <- cellcond$treatment!="-"
dat <- dat[,keep]
cellcond <- cellcond[keep,]

keep <- cellcond$treatment!="empty check"  #not sure what these are. remove
dat <- dat[,keep]
cellcond <- cellcond[keep,]

#one messed up control
# keep <- cellcond$libraryname!="samp_15"
# dat <- dat[,keep]
# cellcond <- cellcond[keep,]




#Number of exonic reads
gene_count <- as.integer(colSums(dat[grep("ENSMUSG",rownames(dat)),]))
cellcond$gene_count <- gene_count 
#Number of detected genes
detected_genes <- as.integer(colSums(dat>0))
#Proportion of mitochondrial reads ------ this one is a bit suspicious...
mt_counts <- as.integer(colSums(na.rm = TRUE,dat[mt_genes$ensembl_gene_id,]))
mt_prop <- mt_counts / gene_count  



findbadcells <- function(forcells=1:ncol(dat),mt_cutoff=0.1){
  # Using scater to assign outliers
  libsize.drop <- isOutlier(gene_count[forcells], nmads=2, type="lower",log=TRUE)
  feature.drop <- isOutlier(detected_genes[forcells], nmads=2, type="lower", log=TRUE)
  mito.drop <- mt_prop[forcells]>mt_cutoff

  todrop <- libsize.drop | feature.drop | mito.drop
  print(dim(todrop))
  totc<-length(todrop)
  print(as.matrix(list(
    total=totc,
    kept=totc-length(which(todrop)),
    drop.by.size=length(which(libsize.drop)),
    drop.by.features=length(which(feature.drop)),
    drop.by.mito=length(which(mito.drop))
  )))

  todroptot <- rep(FALSE, ncol(dat))
  todroptot[(1:ncol(dat))[forcells][todrop]] <- TRUE

  par(mfrow=c(1,2))
  hist(log10(gene_count[forcells]), xlab="Mapped reads (log-scale)", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(gene_count[forcells][libsize.drop])
  abline(v=log10(v), col="red", lty=2, lwd=1.5)
  hist(detected_genes[forcells], xlab="Number of detected genes", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(detected_genes[forcells][feature.drop])
  abline(v=v, col="red", lty=2, lwd=1.5)
  # hist(mt_prop[forcells], xlab="mtDNA %", main="",breaks=50, col="grey80", ylab="Number of cells")
  # abline(v=mt_cutoff, col="red", lty=2, lwd=1.5)
  #dev.off()

  todroptot
}

#Remove crappy libraries
pdf("upstream2/qc_cytoscreen.pdf")
toremove <- findbadcells()
dev.off()
cellcond[toremove,]
keep <- !toremove
dat <- dat[,keep]
cellcond <- cellcond[keep,]




###### Dropped. Seems fairly random, so ok
# well treatment mouse libraryid plate libraryname total_reads gene_count
# 1    A1     empty     1         1     1      samp_1    486.0000        486
# 6    A3     empty     1         3     1      samp_3      0.0000          0
# 8    A5     empty     1         5     1      samp_5    519.0001        519
# 42   D3     mCCL7     1        39     1     samp_39    726.9402        726
# 54   E3    mCCL22     1        51     1     samp_51    792.0000        792


cellcond$filename1 <- sprintf("s%s_S%s_L001_R1_001.fastq.gz",cellcond$libraryid, cellcond$libraryid)
cellcond$filename2 <- sprintf("s%s_S%s_L002_R1_001.fastq.gz",cellcond$libraryid, cellcond$libraryid)

write.csv(cellcond, "cyto2_cellcond.csv")
write.csv(dat,      "cyto2_dat.csv")

# allmd5 <- read.csv("cyto2_dat_md5.csv")
# head(allmd5)
# allmd5 <- allmd5[allmd5$filename %in% c(cellcond$filename1, cellcond$filename2),]
# write.csv(allmd5,"cyto2_dat_md5_some.csv")

######################################################################
############# Prepare for analysis ###################################
######################################################################

ensidFoxp3 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Foxp3"]
ensidRora <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Rora"]


cellcond$depth <- apply(dat,2,sum)
cellcond$treatment <- factor(cellcond$treatment)
cellcond$plate <- factor(cellcond$plate)
cellcond$mouse <- factor(cellcond$mouse)

## Normalize
ncnt <- normalizeDeseqByAll(dat)

av_ncnt <- data.frame(
  ensembl_gene_id=rownames(ncnt),
  mean=apply(ncnt,1,mean))

# mean(ncnt[ensidFoxp3,cellcond$treatment=="wt"])
# mean(ncnt[ensidFoxp3,cellcond$treatment=="our_TGFb"])





######################################################################
############# Cluster conditions #####################################
######################################################################

#Batch-correct normalized counts
batchcorr <- limma::removeBatchEffect(
  #  ncnt,
  log2(1+ncnt),
  batch = cellcond$mouse,
  covariates = cellcond$depth)


keepgenes <- apply(ncnt,1,mean)>0.2
sum(keepgenes)
#nrow(ncnt)

set.seed(0)
ncntred <- batchcorr[keepgenes,]
#ncntred <- dat[keepgenes,]
#nrow(ncntred)
d <- stats::dist(t((ncntred)))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=5, verbose = TRUE, max_iter=5000,dims = 2)
keep <- cellcond$treatment!="foo"
plot(rtsne_out$Y[keep,], pch=16, main='',xlab="",ylab="",cex=0)  
# text(rtsne_out$Y[keep,1], rtsne_out$Y[keep,2],
#      labels=cellcond$treatment[keep],cex=1,   col=tabrora$colplate[keep])
text(rtsne_out$Y[keep,1], rtsne_out$Y[keep,2],
     labels=sprintf("%s_%s",cellcond$treatment[keep],cellcond$libraryid[keep]),
     cex=1,   col=tabrora$colplate[keep])

#wt29, wt_2

######################################################################
############# check Rora #############################################
######################################################################




tabrora <- cellcond
tabrora$rora=ncnt[ensidRora,]

##Color by condition
tabrora$col <- "black"
tabrora$col[tabrora$treatment=="wt"] <- "blue"

##Color by plate
tabrora$colplate <- "black"
tabrora$colplate[tabrora$plate==2] <- "red"
tabrora$colplate[tabrora$plate==4] <- "blue"


tabrora <- tabrora[order(tabrora$rora),]
tabrora <- tabrora[tabrora$mouse==4,]
plot(tabrora$rora, col=tabrora$col)


######################################################################
############# Linear model to compare conditions #####################
######################################################################


########## Fit linear model
datred <- dat[apply(ncnt,1,mean)>0.3,]
nrow(datred)
cellcond$mouse     <- factor(cellcond$mouse)
cellcond$treatment <- factor(cellcond$treatment)
# dds <- DESeqDataSetFromMatrix(countData = round(datred),
#                               colData = cellcond,
#                               design = ~treatment + mouse)
dds <- DESeqDataSetFromMatrix(countData = round(datred),
                              colData = cellcond,
                              design = ~treatment)
dds <- DESeq(dds)


########## Test all conditions vs WT
cyto_de <- list()
for(curcond in unique(cellcond$treatment)){
  if(curcond!="wt"){
    print(curcond)
    v <- results(dds,contrast=c("treatment","wt",curcond))
    v$ensembl_gene_id <- rownames(v)
    #v <- merge(as.data.frame(v),ensconvert)  #better do later!
    v <- as.data.frame(v[order(v$pvalue),])
    cyto_de[[curcond]] <- v
    #v
  }
}


sum(v$pvalue<0.01)  #a measure of phenotype strength
#cyto_de[[curcond]]$

#v[v$mgi_symbol=="Rora",]


sort(table(cellcond$treatment))
#47 controls. all treatments in at least duplicate

########## Assembly FCs and pvalues from the tests
testedcond <- names(cyto_de)
#cyto_de <- list()
v<-cyto_de[[testedcond[1]]]
cyto_fc    <- matrix(ncol=length(testedcond), nrow=nrow(v))
cyto_p     <- matrix(ncol=length(testedcond), nrow=nrow(v))
cyto_fcse  <- matrix(ncol=length(testedcond), nrow=nrow(v))
rownames(cyto_fc)   <- v$ensembl_gene_id
rownames(cyto_p)    <- v$ensembl_gene_id
rownames(cyto_fcse) <- v$ensembl_gene_id
colnames(cyto_fc)   <- testedcond
colnames(cyto_p)    <- testedcond
colnames(cyto_fcse) <- testedcond
#store lfSE too?
for(i in 1:length(testedcond)){
  curcond <- testedcond[i]
  v <- cyto_de[[curcond]]
  cyto_fc  [v$ensembl_gene_id,i] <- v$log2FoldChange
  cyto_p   [v$ensembl_gene_id,i] <- v$pvalue
  cyto_fcse[v$ensembl_gene_id,i] <- v$lfcSE
}

#colnames(cyto_fc) #"empty check" "mSDF1a"      "mCCL7"       "mCCL22"      "hInfb"       "mIL33"       "mTNFa" 

# colnames(cyto_fc) <- c(
#   "mCCL7","mCXCl20","IFNb","mCXCL10","mIL5","mCXCL16","hIL24","mIL4 f row",
#   "mSDF1a","mIL13","mCCL2 rep a","mTNFSF18","mCCL22 rep a","mCXCL9","mTNFa","mIL18",
#   "mIFNa","mIL21","mIL33","mIL10","mCCL2 rep b","mIL4 rep our","mCCL5","mIL6 rep a",
#   "mIFNg","mIL15","mCCL17","mCCL20","mIL1a","rmIL1b","mIL6 rep b","mCCL27",
#   "mCCL25","rmIL2b","mCCL21","mCCL22 rep b","Danazol","mTSLP","mTGFb1","rmTGFb",
#   "mTGFb rep our","mIL2","mIL9","mCCL4","Estrogen_17a")
colnames(cyto_p) <- colnames(cyto_fc)
colnames(cyto_fcse) <- colnames(cyto_fc)

write.csv(cyto_fc,  "out_upstream/cytosceen2_fc.csv")
write.csv(cyto_p,   "out_upstream/cytosceen2_p.csv")
write.csv(cyto_fcse,"out_upstream/cytosceen2_fcse.csv")


cyto_fc   <- read.csv("out_upstream/cytosceen2_fc.csv", check.names = FALSE, stringsAsFactors = FALSE)
cyto_p    <- read.csv("out_upstream/cytosceen2_p.csv", check.names = FALSE, stringsAsFactors = FALSE)
cyto_fcse <- read.csv("out_upstream/cytosceen2_fcse.csv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(cyto_fc)   <- cyto_fc[,1]
rownames(cyto_p)    <- cyto_p[,1]
rownames(cyto_fcse) <- cyto_fcse[,1]
cyto_fc   <- cyto_fc[,-1]
cyto_p    <- cyto_p[,-1]
cyto_fcse <- cyto_fcse[,-1]

######################################################################
############# Plot Rora FC for paper #################################
######################################################################

#verify, what is mouse and human?

### color by family, CCL, IL, etc


### TODO. set FC correct in the input data instead
plotGeneFC <- function(genename, o=NULL){
  
#  genename="Rora"
  geneid <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol==genename]
    
  thefc <- data.frame(
    name=colnames(cyto_fc),
    fc=-as.double(cyto_fc[geneid,]), 
    se=as.double(cyto_fcse[geneid,]),
    stringsAsFactors = FALSE)

  if(is.null(o)){
    o <- order(thefc$fc)
  }
  
  thefc <- thefc[o,]
  thefc$name <- factor(thefc$name,levels = thefc$name)
  thefc$fclow  <- thefc$fc-thefc$se
  thefc$fchigh <- thefc$fc+thefc$se
  
  p<-ggplot(
    thefc,
    aes(x=name, y=fc, ymin=fclow, ymax=fchigh)
  ) + geom_col() + geom_errorbar() + theme_classic() +  
    labs(y=sprintf("%s log2 FC", genename), x=NULL) +
    theme(
      axis.text.x = element_text(angle = (-90), size = 13, face="bold",hjust=0, vjust=1),
      axis.text.y = element_text(angle = (-90), size = 13, face="bold",hjust=0, vjust=1)
    )
  print(p)
  
  o
}




plotGeneFC_2 <- function(genename1, genename2){
  
  #  genename="Rora"
  geneid1 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol==genename1]
  geneid2 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol==genename2]
  
  thefc <- data.frame(
    name=colnames(cyto_fc),
    fc1=-cyto_fc[geneid1,], 
    se1=cyto_fcse[geneid1,],
    fc2=-cyto_fc[geneid2,], 
    se2=cyto_fcse[geneid2,],
    stringsAsFactors = FALSE)
  
  thefc$name <- factor(thefc$name,levels = thefc$name)
  # thefc$fclow  <- thefc$fc-thefc$se
  # thefc$fchigh <- thefc$fc+thefc$se
  
  plot(thefc$fc1, thefc$fc2,cex=0,xlab=genename1, ylab=genename2)
  lines(c(-5,5), c(0,0), col="gray")
  lines(c(0,0), c(-5,5), col="gray")
  text(thefc$fc1, thefc$fc2, labels = thefc$name)
  
}


pdf("out_upstream/rescreen_rora_barplot.pdf")
plotGeneFC("Rora")
dev.off()

#calculate correlation between treatments... actually hclust!

# hdat <- batchcorr
# #hdat <- log(1+dat)
# colnames(hdat) <- sprintf("%s_%s",cellcond$treatment,cellcond$libraryname)
# plot(hclust(dist(t(hdat))))
# 
# 
# #t.test(
#   batchcorr[ensidRora,cellcond$treatment=="wt"]
#   batchcorr[ensidRora,cellcond$treatment=="mCCL22"]
#   #)
# t.test(
#   batchcorr[ensidRora,cellcond$treatment=="wt"],
#   batchcorr[ensidRora,cellcond$treatment=="mTNFa"])
# 
# 
# round(apply(dat,2,sum))
# 
# 
# mouse_gene_len <- read.csv("mouse_gene_length.csv",sep="\t")
# colnames(mouse_gene_len) <- c("ensembl_gene_id","len")
# length(unique(mouse_gene_len$ensembl_gene_id))
# mouse_gene_len <- sqldf("select ensembl_gene_id, avg(len) as len from mouse_gene_len group by ensembl_gene_id")
# rownames(mouse_gene_len) <- mouse_gene_len$ensembl_gene_id
# 
# 
# head(mouse_gene_len)
# ensconvert$ensembl_gene_id
# 
# rpkm <- function(counts, lengths) {
#   rate <- counts / lengths 
#   rate / sum(counts) * 1e6
# }
# 
# calc_tpm <- function(counts, lengths) {
#   rate <- counts / lengths
#   narate <- rate[!is.na(rate)]
#   rate / sum(narate) * 1e6
# }
# 
# avtpm <- calc_tpm(
#   apply(dat,1,sum),
#   mouse_gene_len[rownames(dat),]$len)
# 
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Rora"]]   #14
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Rorc"]]   #9.7
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Gata3"]]  #610
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Foxp3"]]  #2.9
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Tbx21"]]  #63
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Il4"]]    #12
# avtpm[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Il13"]]   #1.9
# 
# 
# 
# 
# 
# 
# ##########
# keepgenes <- which(apply(abs(cyto_fc)-cyto_fcse>0.5,1,sum) > 3)
# red <- cyto_fc[keepgenes,]
# dim(red)
# 
# #can also have an exp level cutoff
# 
# rcor <- cor(red[ensid_rora,], t(red), method = "spearman")
# rcor <- rcor - median(rcor)
# hist(as.double(rcor))
# rcor <- as.data.frame(t(rcor))
# rcor$ensembl_gene_id <- row.names(rcor)
# rcor <- merge(rcor, ensconvert)
# rcor <- rcor[order(rcor$V1),]
# 
# rcor <- rcor[rcor$ensembl_gene_id %in% av_ncnt$ensembl_gene_id[av_ncnt$mean>2],]
# 
# 
# head(rcor,n=50)
# #Stat5b !!!!
# tail(rcor,n=50)
# ## Trak2, B3gnt3  Bach1
# ## il6ra !! Bach2  
# #############
# 
# 
# 
# ensconvert$ensembl_gene_id
# rcor
# #rcor
# #rcor$foo <- 1
# names(rcor[rcor < -0.5])
# #rcor <- rcor[order(rcor)]  #sort(rcor)
# head(rcor)
# #rcor <- rcor[abs(rcor)>]
# sort(rcor)
# 
# head(red)
# 
# 
# pdf("out_upstream/keep_cytokine_rora_barplot.pdf")
# orderrora <- plotGeneFC("Rora")
# dev.off()
# pdf("out_upstream/cytokine_foxp3_barplot_ordered_by_rora.pdf")
# plotGeneFC("Foxp3",o = orderrora) 
# dev.off()
# pdf("out_upstream/keep_cytokine_foxp3_barplot.pdf")
# plotGeneFC("Foxp3") 
# dev.off()
# pdf("out_upstream/keep_cytokine_tbx21_barplot.pdf")
# plotGeneFC("Tbx21") 
# dev.off()
# pdf("out_upstream/keep_cytokine_gata3_barplot.pdf")
# plotGeneFC("Gata3") 
# dev.off()
# pdf("out_upstream/keep_cytokine_rorc_barplot.pdf")
# plotGeneFC("Rorc") 
# dev.off()
# 
# 
# plotGeneFC("Il6ra") 
# 
# 
# 
# 
# pdf("out_upstream/cytokine_bach2_barplot.pdf")
# plotGeneFC("Bach2")
# dev.off()
# 
# pdf("out_upstream/cytokine_tnfrsf25_barplot.pdf")
# plotGeneFC("Tnfrsf25")
# dev.off()
# 
# 
# ######## for bm/jp
# pdf("out_upstream/cytokine_xbp1_barplot.pdf")
# plotGeneFC("Xbp1")
# dev.off()
# pdf("out_upstream/cytokine_ern1_barplot.pdf")
# plotGeneFC("Ern1")
# dev.off()
# pdf("out_upstream/cytokine_cyp11a1_barplot.pdf")
# plotGeneFC("Cyp11a1")
# dev.off()
# 
# 
# 
# pdf("out_upstream/cytokine_St6galnac3_barplot.pdf")
# plotGeneFC("St6galnac3")
# #plotGeneFC("St6galnac3",o = orderrora)
# dev.off()
# 
# 
# pdf("out_upstream/cytokine_Cxcr6_barplot.pdf")
# plotGeneFC("Cxcr6")
# dev.off()
# 
# 
# pdf("out_upstream/cytokine_Nebl_barplot.pdf")
# plotGeneFC("Nebl")
# dev.off()
# 
# 
# 
# # pdf("out_upstream/circ_cytokine_arntl.pdf")
# # plotGeneFC("Arntl")  
# # dev.off()
# 
# pdf("out_upstream/circ_cytokine_arntl.pdf")
# plotGeneFC("Arntl")  
# dev.off()
# 
# pdf("out_upstream/circ_cytokine_ccr7.pdf")
# plotGeneFC("Ccr7")  
# dev.off()
# 
# 
# plotGeneFC("S1pr1")
# 
# plotGeneFC("Per1")    #IL21, IL9, Tnfrsf18
# plotGeneFC("Per2")    #tricky
# plotGeneFC("Per3")    #IL5, IL9, CCL2
# plotGeneFC("Clock")   #tricky
# plotGeneFC("Nr1d1")   #crappy intervals
# plotGeneFC("Nr1d2")   #crappy intervals
# plotGeneFC("Creb1")   #-Ccl22, maybe more
# plotGeneFC("Mtor")
# 
# 
# 
# 
# 
# ### should do a GO analysis of each perturbation
# 
# plotGeneFC("Calcrl")  
# plotGeneFC("Calca")  
# 
# 
# 
# 
# pdf("out_upstream/cytokine_Prdm1_barplot.pdf")
# plotGeneFC("Prdm1")
# dev.off()
# 
# 
# 
# ######### Overlap with DE axis Th1/Th2/Treg. where do they map?
# 
# 
# plotGeneFC("",orderrora)
# 
# 
# ### we should make a website where you can check your fawourite gene





######################################################################
### Compare the treatments by correlation ############################
######################################################################

# 
# #should only look at genes which change properly
# 
# 
# ### list of DE genes, Rora overexp?
# 
# rownames(de_overexp) <- de_overexp$ensembl_gene_id
# cyto_fc_rora <- cbind(cyto_fc, de_overexp[rownames(cyto_fc),]$fc_overexp)
# cyto_fc_rora <- cyto_fc_rora[apply(!is.na(cyto_fc_rora),1,all),]
# 
# #cyto_fc_de <- cyto_fc_rora
# cyto_fc_de <- cyto_fc[,!(colnames(cyto_fc) %in% c("mIL6 rep b","mIL21","Estrogen_17a","mCCL4"))]
#   #cyto_fc[apply(cyto_fc>0.5,1,any),]
# 
# ##Only keep the most relevant genes
# cyto_fc_de <- cyto_fc_de[apply(abs(cyto_fc_de)>1,1,any),]
# dim(cyto_fc_de)
# 
# ##Get rid of the diagonal
# tcorr <- cor(cyto_fc_de)
# for(i in 1:nrow(tcorr)){
#   tcorr[i,i] <- 0
# }
# ##Get rid of technical bias toward positive correlation
# tcorr <- cor(tcorr)
# 
# ##Normalize correlations to get "strength". actually, remove average?
# sort(apply(tcorr,1,mean))
# 
# ensidIL6 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Il6"]
# ensidIL4 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Il4"]
# 
# ensidXXX <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Il17a"]
# ensidXXX <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Gata3"]
# 
# 
# ### Put rora FC on the side!!
# sidecol_data <- cyto_fc_de[ensidRora,]
# #sidecol_data <- cyto_fc_de[ensidIL6,]
# sidecol_data <- cyto_fc_de[ensidXXX,]
# 
# 
# sidecol_index <- as.integer(sidecol_data/max(abs(sidecol_data))*500)+500
# #sidecol_index <- as.integer(sidecol_data/max(abs(sidecol_data))*1000)+500
# sidecol_index[sidecol_index<1]    <- 1
# sidecol_index[sidecol_index>1000] <- 1000
# sidecol <- colorRampPalette(c("green", "black", "red"))(n = 1000)[sidecol_index]
# 
# #pdf("out_upstream/corcor_treatments.pdf")
# plot.new()
# heatmap.2(tcorr*5,
#           ColSideColors = sidecol,
#           RowSideColors = sidecol,
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           col=colorRampPalette(c("green", "black", "red"))(n = 1000)
#           # margins =c(12,9),     # widens margins around plot
# )            # turn off column clustering
# #dev.off()
# 
# 
# ### Which genes are commonly increasing in upregulating Rora?
# ### Which genes are commonly decreasing in upregulating Rora?
# 
# #plus: mCCL5, mIL15, mIL33, mCCL7
# #negative: mCCL22 rep a, mSDF1a, IFNb, mCCL25
# 
# correlating <- cbind(
#   cyto_fc[,c("mCCL5", "mIL15", "mIL33", "mCCL7")],
#   -cyto_fc[,c("mCCL22 rep a", "mSDF1a", "IFNb", "mCCL25")])
# 
# correlating <- cbind(
#   cyto_fc[,c("mCCL5", "mIL15", "mIL33")],  #, "mCCL7"
#   -cyto_fc[,c("mCCL22 rep a", "mSDF1a", "IFNb")])  #, "mCCL25"
# 
# #av_corr <- sort(apply(correlating,1,mean))
# 
# av_corr <- merge(
#   ensconvert,  
#   data.frame(
#     ensembl_gene_id=names(av_corr),
#     m=av_corr))
# 
# av_corr <- av_corr[order(av_corr$m),]
# av_corr[1:100,]
# ### Rora    Il4i1  Tnfrsf21    st3gal6     Ppargc1a    
# 
# 
# av_corr <- av_corr[order(av_corr$m, decreasing = TRUE),]
# av_corr[1:1000,]
# ### mt-Apt6   il11ra1   
# 
# ### could do a GO-analysis here
# 
# ### not just FC - should focus on good p-value genes
# 

######################################################################
### Fancy way of comparing cytokine effects, using marker genes ######
######################################################################
# 
# ### Get marker genes from Th-express. 
# d <-         read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_naive.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE]
# d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_Th1.txt",  sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
# d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_Th17.txt", sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
# d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_iTreg.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
# d <- cbind(d,read.csv("/home/mahogny/Dropbox/applyPI/crpaper/code/analysis/thexpress/Th2_vs_nTreg.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
# thep <- d[,c(2,4,6,8,10)]
# thefc <- d[,c(1,3,5,7,9)]
# colnames(thep)<-c("Naive/Th0","Th1","Th17","iTreg","nTreg")
# colnames(thefc)<-c("Naive/Th0","Th1","Th17","iTreg","nTreg")
# 
# for(i in 1:ncol(thep)){
#   thep[is.na(thep[,i]),i] <- 1
# }
# for(i in 1:ncol(thefc)){
#   thefc[is.na(thefc[,i]),i] <- 0
# }
# thmarkergenes <- 
#   intersect(
#     names(which(apply(thep<1e-17 & abs(thefc)>5,1,any))),
#     rownames(cyto_fc)
#     )
# #thmarkergenes_sym <- sort(togenesym(thmarkergenes))
# 
# 
# #thmarkergenes %in% rownames(cyto_fc)
# 
# 
# cyto_fc_marker <- cyto_fc[thmarkergenes,]
# thefc_marker <- thefc[thmarkergenes,]
# 
# 
# #cor(cyto_fc)
# 
# 
# 
# 
# 
# v <- cor(thefc_marker$iTreg, cyto_fc_marker, method = "spearman")[1,]  #TGFb at top
# v <- v[order(v)]
# v
# 
# v <- cor(thefc_marker$nTreg, cyto_fc_marker, method = "spearman")[1,]   #nTreg... estrogen at top, interesting!
# v <- v[order(v)]
# v
# 
# v <- cor(thefc_marker$Th17, cyto_fc_marker, method = "spearman")[1,]  
# v <- v[order(v)]
# v
# 
# v <- cor(thefc_marker$Th1, cyto_fc_marker, method = "spearman")[1,]  
# v <- v[order(v)]
# v
# 
# v <- cor(thefc_marker$`Naive/Th0`, cyto_fc_marker)[1,]  
# v <- v[order(v)]
# v
# 
# 
# 
# 
# thefc$`Naive/Th0`
# 
# thmarkergenes
# 
# thmarkergenes_sym
# 
# 
# 
# 
# #thmarkergenes
# # c("Foxp3","Tbx21","Gata3") %in% thmarkergenes_sym
# # length(thmarkergenes_sym)
# 

######################################################################
### What do the genes in the groups have in common? ##################
######################################################################


# groupA <- c(
#   "mTNFSF18",
#   "mIL1a",
#   "mTSLP",
#   "IFNb",
#   "rmIL1b",
#   "mIL4 rep our",
#   "mTGFb1",
#   "mCCL21",
#   "mIL6 rep a",
#   "mIL10",
#   "mIFNa",
#   "rmTGFb",
#   "mCCL22 rep a",
#   "mCCL17",
#   "mTGFb rep our",
#   "mIFNg",
#   "mIL33",
#   "mCCL20",
#   "rmIL2b"
#   )
# 
# groupB <- c(
#   "Danazol",
#   "mSDF1a",
#   "mCXCL16",
#   "mIL13",
#   "mCCL7",
#   "mCCL22 rep b",
#   "mCCL2 rep a",
#   "mIL5",
#   "mCCL27",
#   "mCXCl20",
#   "mIL4 f row",
#   "mIL2",
#   "mCCL5"
# )
# 
# 
# 
# 
# groupA %in% colnames(cyto_fc)
# groupB %in% colnames(cyto_fc)
# 
# diff <- apply(cyto_fc[,groupB],1,mean)-apply(cyto_fc[,groupA],1,mean)
# 
# ##########
# ########## Compare A vs B
# ##########
# diff <- apply(cyto_fc[,groupB],1,mean)-apply(cyto_fc[,groupA],1,mean)
# diff <- sort(abs(diff), decreasing = TRUE)
# 
# stopgosym(
#   togenesym(names(diff[1:500])),
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% names(diff)])
#   )
# 
# 
# ##########
# ########## Compare A vs WT
# ##########
# diff <- apply(cyto_fc[,groupA],1,mean)
# diff <- sort(abs(diff), decreasing = TRUE)
# 
# stopgosym(
#   togenesym(names(diff[1:300])),
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% names(diff)])
# )
# #response to cytokine       584         580   575.68  0.0771
# 
# 
# 
# 
# 
# ##########
# ########## Compare B vs WT
# ##########
# diff <- apply(cyto_fc[,groupA],1,mean)
# diff <- sort(abs(diff), decreasing = TRUE)
# 
# stopgosym(
#   togenesym(names(diff[1:300])),
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% names(diff)])
# )
# 
# 
# 
# 
# ######################################################################
# ### Full DE test of group A ##########################################
# ######################################################################
# 
# groupAold <- c(
#   "mtnfsf18",
#   "mil1a",
#   "mtslp",
#   "ifnb",
#   "rmIL1b",
#   "_il4", #maybe #"mIL4 rep our",
#   "mouse_TGFb1",   #maybe mTGFb
#   "21", ##  ccl21
#   "mIL6",
#   "mIL10",
#   "ifna",
#   "rmTGFb",
#   "22",
#   "17",
#   "our_TGFb",
#   "IFNg",
#   "mil33",
#   "20",
#   "rmIL2b"
# )
# 
# cellcond_a <- cellcond
# cellcond_a$group <- ""
# cellcond_a$group[cellcond_a$treatment %in% groupAold] <- "t"
# cellcond_a$group[cellcond_a$treatment %in% "wt"] <- "c"
# cellcond_a$group <- factor(cellcond_a$group)
# sort(cellcond_a$treatment)
# 
# keep <- which(cellcond_a$treatment!="")
# 
# 
# ########## Fit linear model
# datred <- dat[apply(ncnt,1,mean)>0.1,keep]
# nrow(datred)
# dds <- DESeqDataSetFromMatrix(countData = round(datred),
#                               colData = cellcond_a[keep,],
#                               design = ~group + mouse)
# dds <- DESeq(dds)
# de_groupa <- results(dds,contrast=c("group","c","t"))
# de_groupa$ensembl_gene_id <- rownames(de_groupa)
# de_groupa <- as.data.frame(de_groupa[order(de_groupa$pvalue),])
# 
# 
# # stopgo(
# #   de_groupa$ensembl_gene_id[1:100],
# #   de_groupa$ensembl_gene_id, 
# #   ID="EnsemblID")
# 
# stopgosym(
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_groupa$ensembl_gene_id[1:100]]),
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_groupa$ensembl_gene_id]))
# 
# 
# 
# ####################
# #################
# 
# 
# groupAsym <- c(
#   "Tnfsf18",
#   "Il1a",
#   "Tslp",
#   "Ifnb",
#   "Il1b",
#   "Il4",
#   "Tgfb1",
#   "Ccl21",
#   "Il6",
#   "Il10",
#   "Ifna",
#   "Tgfb",
#   "Ccl22",
#   "Ccl17",
#   "Tgfb",
#   "Ifng",
#   "Il33",
#   "Ccl20",
#   "Il2b"
# )
# 
# groupBsym <- c(
#   "Sdf1a",
#   "Cxcl16",
#   "Il13",
#   "Ccl7",
#   "Ccl22",
#   "Ccl2",
#   "Il5",
#   "Ccl27",
#   "Cxcl20",
#   "Il4",
#   "Il2",
#   "Ccl5"
# )
# 
# groupABsym <- unique(c(groupAsym, groupBsym))
# 
# stopgosym(groupAsym, groupABsym)
# stopgosym(groupBsym, groupABsym)   #This one pulls out "regulation of interleukin-6 production", p=0.0077
# 
# 
# 
# 
# 
# 
# ######################################################################
# ### Full DE test of group A ##########################################
# ######################################################################
# 
# 
# groupBold <- c(
#   "danazol",
#   "sdf1a",
#   "mcxcl16",
#   "il13",
#   "mCcl7",
#   "22",
#   "2",
#   "IL5",
#   "27",
#   "mcxcl20",
#   "our_IL2",
#   "5"
# )
# 
# sort(unique(cellcond$treatment))
# 
# 
# cellcond_b <- cellcond
# cellcond_b$group <- ""
# cellcond_b$group[cellcond_b$treatment %in% groupBold] <- "t"
# cellcond_b$group[cellcond_b$treatment %in% "wt"] <- "c"
# cellcond_b$group <- factor(cellcond_b$group)
# #sort(cellcond_b$treatment)
# 
# keep <- which(cellcond_b$treatment!="")
# 
# 
# ########## Fit linear model
# datred <- dat[apply(ncnt,1,mean)>0.1,keep]
# nrow(datred)
# dds <- DESeqDataSetFromMatrix(countData = round(datred),
#                               colData = cellcond_a[keep,],
#                               design = ~group + mouse)
# dds <- DESeq(dds)
# de_groupb <- results(dds,contrast=c("group","c","t"))
# de_groupb$ensembl_gene_id <- rownames(de_groupb)
# de_groupb <- as.data.frame(de_groupb[order(de_groupb$pvalue),])
# 
# 
# stopgosym(
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_groupb$ensembl_gene_id[1:100]]),
#   unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_groupb$ensembl_gene_id]))
# 
# 
# 
# 
# # Estrogen suppresses T cell TNF production by regulating T cell differentiation and activity in the bone marrow,
# # thymus, and peripheral lymphoid organs. In the bone marrow, estrogen downregulates the proliferation of hematopoietic
# # stem cells through an IL-7 dependent mechanism.[18]
# 




#########################################################################
###################### Compare with overexp #############################
#########################################################################


#### compare p value
overlap_overexp_cyto <- list()
for(i in 1:ncol(cyto_p)){
  overlap_overexp_cyto[[colnames(cyto_p)[i]]] <- intersect(
    rownames(cyto_p)[order(cyto_p[,i])][1:200],
    de_merge$ensembl_gene_id[order(de_merge$p_overexp)][1:300]    #p_th17
  )
}
overlap_overexp_cyto
sort(as.data.frame(lapply(overlap_overexp_cyto, length)))


#### compare absolute FC
overlap_overexp_cyto <- list()
for(i in 1:ncol(cyto_p)){
  overlap_overexp_cyto[[colnames(cyto_fc)[i]]] <- intersect(
    rownames(cyto_fc)[order(abs(cyto_fc[,i]), decreasing = TRUE)][1:200],
    de_merge$ensembl_gene_id[order(abs(de_merge$fc_th17), decreasing = TRUE)][1:200]  
  )
}
#overlap_overexp_cyto

sort(as.data.frame(lapply(overlap_overexp_cyto, length)))

cyto_fc$mCCL7

out <- c()
for(i in 1:ncol(cyto_fc)){
  x<-merge(
    data.frame(
      ensembl_gene_id=rownames(cyto_fc),
      fc_cyto=cyto_fc[,i]
    ),
    data.frame(
      ensembl_gene_id=de_merge$ensembl_gene_id,
      fc_merge=de_merge$fc_overexp  
    )
  )
  x <- x[abs(x$fc_merge)>2,]   #2 for out_2
  out <- c(out,cor((x$fc_cyto),(x$fc_merge)))
}
names(out) <- colnames(cyto_fc)
out
out <- sort(out)
#mean(sign(x$fc_cyto)*sign(x$fc_merge))



#out_4 <- sort(out)
#out_3 <- sort(out)
#out_2 <- sort(out)
## cytokine inducer should have positive correlation with overexpression of rora, to be thought to induce rora

nout <- data.frame(
  name=names(out_2),
  cor=out_2,
  stringsAsFactors = FALSE
  )

nout$name <- factor(nout$name,levels = nout$name)

p<-ggplot(
  nout,
  aes(x=name, y=cor)
) + geom_col() + theme_classic() +  
  labs(y=sprintf("%s cor FC", nout$name), x=NULL) +
  theme(
    axis.text.x = element_text(angle = (-90), size = 13, face="bold",hjust=0, vjust=1),
    axis.text.y = element_text(angle = (-90), size = 13, face="bold",hjust=0, vjust=1)
  )
pdf("cytocor_repeat.pdf")
print(p)
dev.off()


#########################################################################
###################### Compare with overexp, FC #############################
#########################################################################



#### compare p value
overlap_overexp_cyto <- list()
for(i in 1:ncol(cyto_p)){
  overlap_overexp_cyto[[colnames(cyto_p)[i]]] <- intersect(
    rownames(cyto_p)[order(cyto_p[,i])][1:200],
    de_merge$ensembl_gene_id[order(de_merge$p_overexp)][1:300]    #p_th17
  )
}
overlap_overexp_cyto
sort(as.data.frame(lapply(overlap_overexp_cyto, length)))


#### compare absolute FC
overlap_overexp_cyto <- list()
for(i in 1:ncol(cyto_p)){
  overlap_overexp_cyto[[colnames(cyto_fc)[i]]] <- intersect(
    rownames(cyto_fc)[order(abs(cyto_fc[,i]), decreasing = TRUE)][1:200],
    de_merge$ensembl_gene_id[order(abs(de_merge$fc_th17), decreasing = TRUE)][1:200]  
  )
}
#overlap_overexp_cyto

sort(as.data.frame(lapply(overlap_overexp_cyto, length)))

out <- c()
for(i in 1:ncol(cyto_fc)){
  x<-merge(
    data.frame(
      ensembl_gene_id=rownames(cyto_fc),
      fc_cyto=cyto_fc[,i]
    ),
    data.frame(
      ensembl_gene_id=de_merge$ensembl_gene_id,
      fc_merge=de_merge$fc_overexp  
    )
  )
  x <- x[abs(x$fc_merge)>2,]     # used 2!!!
  foo <- cor.test(x$fc_cyto,x$fc_merge)
  out <- rbind(out,c(foo$estimate, foo$conf.int))
  # foo <- cor(x$fc_cyto,x$fc_merge)
  # out <- rbind(out,c(foo, 0,0))
}
rownames(out) <- colnames(cyto_fc)
colnames(out) <- c("est","lower","upper")
out <- as.data.frame(out)
out
out <- out[order(out$est),]  #sort(out)
#mean(sign(x$fc_cyto)*sign(x$fc_merge))

## cytokine inducer should have positive correlation with overexpression of rora, to be thought to induce rora

nout <- data.frame(
  name=rownames(out),
  cor=out$est,
  lower=out$lower,
  upper=out$upper,
  stringsAsFactors = FALSE
)

nout$name <- factor(nout$name,levels = nout$name)

p<-ggplot(
  nout,
  aes(x=name, y=cor, ymin=lower, ymax=upper)
) + geom_col() + geom_errorbar() + theme_classic() +  
  #labs(y=sprintf("%s cor FC", nout$name), x=NULL) +
  theme(
    axis.text.x = element_text(angle = (-90), size = 13, face="bold",hjust=0, vjust=1),
    axis.text.y = element_text(angle = (-90), size = 13, face="bold",hjust=0, vjust=1)
  )
print(p)
pdf("cytocor_repeat.pdf")
print(p)
dev.off()
