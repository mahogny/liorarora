library("DESeq2")
library("sqldf")
library("matrixStats")
library("statmod")
library("stringr")
library("gplots")
library("limma")
library(biomaRt)
library(org.Mm.eg.db)
library(scater)

##########################################################################
############### Common functions #########################################
##########################################################################

###############################
## Normalization to take genes to same "level"
normalizegenes <- function(genenorm){
  for(i in 1:nrow(genenorm)){
    thesum <- mean(genenorm[i,])
    if(thesum!=0)
      #genenorm[i,] <- genenorm[i,]/thesum
      genenorm[i,] <- genenorm[i,]-thesum
  }
  return(genenorm)
}


normalizegenesSD <- function(genenorm){
  for(i in 1:nrow(genenorm)){
    thesum <- sd(genenorm[i,])
    if(thesum!=0)
      genenorm[i,] <- genenorm[i,]/thesum
  }
  return(genenorm)
}






###############################
## replace row name gene name, with symbolic gene name
ensembl2gene <- read.csv("ensembl2genename.txt",stringsAsFactors=FALSE)
wrealgenename <- function(dat) {
  dat2 <- data.frame(Ensembl.Gene.ID=rownames(dat),stringsAsFactors = FALSE)
  x<-merge(dat2, ensembl2gene, all.x=TRUE)
  thena<-is.na(x[,2])
  x[thena,2] <- x[thena,1]
  rownames(dat)<-x[,2]
  dat
}



###############################
## Exclude genes having no variance. Convenient for working with large genesets and cor() etc
rednovar <- function(y){
  keep <- logical(nrow(y))
  for(i in 1:nrow(y)){
    keep[i] = sd(y[i,])!=0
  }  
  y<-y[which(keep),]
  return(y)
}


###############################
## Plot correlation between genes, but with those having no variation removed
## Diagonal set to 0 for visibility
corgenered <- function(y){
  y<-rednovar(y)
  thecor <- cor(t(y),method="spearman")
  for(i in 1:ncol(thecor))
    thecor[i,i]<-0
  heatmap(thecor)
  return()
}



######################################################################
### Read ensembl #####################################################
######################################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]

#Convert to gene symbols, or retain ID if no symbol
ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart=ensembl),stringsAsFactors=FALSE)
togenesym <- function(ensid){
  out <- c()#terrible speed...
  for(n in ensid){
    ids <- which(ensconvert$ensembl_gene_id==n)
    if(length(ids)==0)
      more <- n
    else
      more <- ensconvert$mgi_symbol[ids]
    out <- c(out, more)    
  }
  #ensconvert$mgi_symbol[which(ensconvert$ensembl_gene_id %in% ensid)]
  out
}

## Convert to ensembl IDs
toensid <- function(geneid){
  out <- c()#terrible speed...
  for(n in geneid){
    out <- c(out, ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol==n)])    
  }
  out
}

ensidPtma <- toensid("Ptma")#ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol=="Ptma")]

nametogenesym<-function(x){
  names(x)<-togenesym(names(x))
  x
}

mtogenesym<-function(x){
  colnames(x)<-togenesym(colnames(x))
  rownames(x)<-togenesym(rownames(x))
  x
}



##########################################################################
################### Read data ############################################
##########################################################################


###############################
## Read data from a single count file
readrawcount <- function(countfile, mappingfile){
  wellmapping <- read.csv(mappingfile, header=FALSE,stringsAsFactors=FALSE)[,1]
  wellmapping <- unlist(lapply(wellmapping,function(x) substr(x,20,str_locate(x,"tagged")[1]-3)))
  
  cnt_schisto <- read.table(countfile,sep="\t",header=TRUE)
  rownames(cnt_schisto) <- cnt_schisto[,1]
  cnt_schisto<-cnt_schisto[,-1]
  cnt_schisto_index<-as.double(unlist(lapply(colnames(cnt_schisto),function(x) substr(x,24,str_locate(x,"counts")[1]-2))))
  # cnt_schisto <- cnt_schisto[-grep("_",rownames(cnt_schisto)),] #remove rows with __
  #"__no_feature"     "__ambiguous"   "__too_low_aQual"  "__not_aligned"  "__alignment_not_unique"
  colnames(cnt_schisto) <- wellmapping[cnt_schisto_index]  
  
  return(cnt_schisto)
}


###############################
## Read data from a directory of count files, one per cell
readrawcountmulti <- function(countfile, mappingfile){
  #  countfile="set2"
  #  mappingfile="mappingSet2.csv"
  wellmapping <- read.csv(mappingfile, header=FALSE,stringsAsFactors=FALSE)[,1]
  wellmapping <- unlist(lapply(wellmapping,function(x) substr(x,20,str_locate(x,"tagged")[1]-3)))
  
  flist<-list.files(countfile)
  cnt_schisto2 <- read.table(sprintf("%s/%s",countfile,flist[1]))
  row.names(cnt_schisto2) <- cnt_schisto2[,1]
  cnt_schisto2<-cnt_schisto2[,c(),drop=FALSE]
  for(fname in flist){
    onef <-sprintf("%s/%s",countfile,fname)
    print(onef)
    cnt_schisto2 <- cbind(cnt_schisto2, read.table(onef)[,2,drop=FALSE])
  }
  #cnt_schisto2 <- cnt_schisto2[-grep("_",rownames(cnt_schisto2)),] #remove rows with __
  colnames(cnt_schisto2)<-flist
  
  cnt_schisto_index<-as.double(unlist(lapply(colnames(cnt_schisto2),function(x) substr(x,1,str_locate(x,"counts")[1]-2))))
  
  #cnt_schisto <- cnt_schisto2[-grep("_",rownames(cnt_schisto2)),] #remove rows with __
  #"__no_feature"     "__ambiguous"   "__too_low_aQual"  "__not_aligned"  "__alignment_not_unique"
  cnt_schisto <- cnt_schisto2
  colnames(cnt_schisto) <- wellmapping[cnt_schisto_index]  
  
  return(cnt_schisto)
}




cnt_file1 <- readrawcountmulti("schisto_c1/set1","schisto_c1/mappingSet1.csv")
cnt_file2 <- readrawcountmulti("schisto_c1/set2","schisto_c1/mappingSet2.csv")

## Concatenate datasets
cnt_schisto <- cbind(cnt_file1,cnt_file2)
#all(rownames(cnt_file1)==rownames(cnt_file2))  #Sanity check
#ncol(cnt_schisto)

### set up a matrix describing the data
cellcondition <- data.frame(isgood=rep(TRUE, ncol(cnt_schisto)))
cellcondition$isMLN <- FALSE
cellcondition$isMLN[grep("MLN",colnames(cnt_schisto))] <- TRUE
cellcondition$isSPL <- FALSE
cellcondition$isSPL[grep("SPL",colnames(cnt_schisto))] <- TRUE

## Ensure data is in format of a double matrix
cnt_schisto <- as.matrix(cnt_schisto)
class(cnt_schisto)<-"double"


##########################################################################
################### RNAseq QC and normalization ##########################
##########################################################################

dat <- cnt_schisto

#Number of exonic reads
gene_count <- colSums(dat[grep("ENSMUSG",rownames(dat)),])
#Number of detected genes
detected_genes <- colSums(dat>0)
#Proportion of mitochondrial reads
mt_counts <- colSums(dat[mt_genes$ensembl_gene_id,])
mt_prop <- mt_counts / gene_count  

findbadcells <- function(forcells,mt_cutoff=0.1){
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
  
  par(mfrow=c(1,3))
  hist(log10(gene_count[forcells]), xlab="Mapped reads (log-scale)", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(gene_count[forcells][libsize.drop])
  abline(v=log10(v), col="red", lty=2, lwd=1.5) 
  hist(detected_genes[forcells], xlab="Number of detected genes", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(detected_genes[forcells][feature.drop])
  abline(v=v, col="red", lty=2, lwd=1.5) 
  hist(mt_prop[forcells], xlab="mtDNA %", main="",breaks=50, col="grey80", ylab="Number of cells")
  abline(v=mt_cutoff, col="red", lty=2, lwd=1.5)
  #dev.off()
  
  todroptot
}

png("out.schisto/QC.png",width=800)
x<-findbadcells(1:nrow(dat))
#as.double(gene_count[x]) #read counts of dropped genes
cellcondition$isgood[x]<-FALSE
dev.off()

### Remove cells with too little coverage as these will just disturb later statistics
cellcondition$isgood[apply(cnt_schisto,2,sum)<1.5e6] <- FALSE


###############################
## Normalization by deseq size factors
normalizeDeseqByGenes <- function(cnt) t(t(cnt)/estimateSizeFactorsForMatrix(cnt[-grep("ERC",rownames(cnt)),]))
normalizeDeseqByERC   <- function(cnt) t(t(cnt)/estimateSizeFactorsForMatrix(cnt[grep("ERC",rownames(cnt)),]))


## Normalize
ncnt_schisto <- normalizeDeseqByGenes(cnt_schisto)

## Store some intermediate tables
# write.table(cnt_schisto, "liora_20140604_schisto.123.txt")
# write.table(ncnt_schisto, "liora_20140604_schisto.123.normalized.txt")

#Translate to common names
tncnt_schisto <- wrealgenename(ncnt_schisto)

#write.table(cnt_schisto, "liora_schisto_raw_data.csv",sep = ",", qmethod = "double")



##########################################################################
################### Define gene expression cut-offs ######################
##########################################################################


library(mixtools)

gene_cutoff<-list()
# cutoffsfor <- rownames(ncnt_schisto)[which(rownames(ncnt_schisto) %in% toensid(c("Rora","Gata3")))]
# for(c in rownames(ncnt_schisto)){
#   print(c)
#   #mean(normalmixEM(log10(as.double(tblsorted[,"RORA"])))$mu)
#   trans <- log10(1+as.double(ncnt_schisto[c,]))
#   table(trans)
#   if(sd(trans)==0){
#     gene_cutoff[[c]] <- 0
#   } else {
#     gene_cutoff[[c]] <- 0
#     try({
#       gene_cutoff[[c]] <- mean(normalmixEM(trans,verb = FALSE)$mu)
#     })
#   }
# }

hist(log10(1+as.double(ncnt_schisto[toensid("Cd4"),])),breaks=20)

#after log10(1+...)
gene_cutoff[["Rora"]] <- 0.5
gene_cutoff[["Sell"]] <- 0.5
gene_cutoff[["Rora"]] <- 0.5
gene_cutoff[["Cd4"]]  <- 0.5
gene_cutoff[["Gata3"]]  <- 0.5
gene_cutoff[["Foxp3"]]  <- 0.5
gene_cutoff[["Tbx21"]]  <- 0.5

#roughly 3 normalized reads

cellcondition$highRora  <- log10(1+as.double(ncnt_schisto[toensid("Rora"),])) > gene_cutoff[["Rora"]]
cellcondition$highGata3 <- log10(1+as.double(ncnt_schisto[toensid("Gata3"),])) > gene_cutoff[["Gata3"]]
cellcondition$highFoxp3 <- log10(1+as.double(ncnt_schisto[toensid("Foxp3"),])) > gene_cutoff[["Foxp3"]]
cellcondition$highTbx21 <- log10(1+as.double(ncnt_schisto[toensid("Tbx21"),])) > gene_cutoff[["Tbx21"]]

cellcondition$highCytokine <- 
  log10(1+as.double(ncnt_schisto[toensid("Il4"),])) > 0.5 |
  log10(1+as.double(ncnt_schisto[toensid("Il10"),])) > 0.5 |
  log10(1+as.double(ncnt_schisto[toensid("Il13"),])) > 0.5

cellcondition$isCD4 <- 
  log10(1+as.double(ncnt_schisto[toensid("Cd4"),])) > 0.5 | 
  log10(1+as.double(ncnt_schisto[toensid("Cd3e"),])) > 0.5



##########################################################################
################### Find genes with real variation #######################
##########################################################################

###############################
## Estimate technical noise
pickAboveTechNoise <- function(cnt_schisto, 
                               pvalue=0.0002){
  sf <- estimateSizeFactorsForMatrix(cnt_schisto)
  ncnt_schisto <- t(t(cnt_schisto)/estimateSizeFactorsForMatrix(cnt_schisto))
  
  meansHeLa <- rowMeans( ncnt_schisto )
  varsHeLa <- rowVars( ncnt_schisto )
  maxHeLa <- rowMaxs( ncnt_schisto )
  cv2HeLa <- varsHeLa / meansHeLa^2
  minMeanForFit <- unname( quantile( meansHeLa[ which( cv2HeLa > .3 ) ], .95 ) )  
  #plot(log(ncnt_schisto[,1]),log(ncnt_schisto[,4]))
  #plot(log(cnt_schisto2[,1]),log(cnt_schisto2[,4]))
  
  useForFit <- meansHeLa >= minMeanForFit
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansHeLa[useForFit] ),
                     cv2HeLa[useForFit] )
  fit$coefficients
  xi <- mean( 1 / sf )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  
  dof <- ncol(cnt_schisto) - 1
  
  ids<-which(cv2HeLa > ((xi+a1)/meansHeLa + a0 )* qchisq(1-pvalue, dof)/dof)
  
  idserc <- grep("ERC",rownames(cnt_schisto))
  
  # Prepare the plot (scales, grid, labels, etc.)
  plot( meansHeLa, cv2HeLa, pch=20, cex=.2, col="gray", xaxt="n", yaxt="n",
        log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .05, 100 ),
        xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                         expression(10^4), expression(10^5) ) )
  axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0" )
  # Add the data points
  #points( meansHeLa, cv2HeLa, pch=20, cex=.2, col="gray")
  points( meansHeLa[ids], cv2HeLa[ids], pch=20, cex=.2, col="green")
  points( meansHeLa[idserc], cv2HeLa[idserc], pch=20, cex=.8, col="red")
  # Plot the fitted curve
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
  # Plot quantile lines around the fit
  lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, dof ) / dof, col="#FF000080", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, dof ) / dof, col="#FF000080", lwd=2, lty="dashed" )
  return(ids)
}


## Pick relevant genes
notnoisy <- pickAboveTechNoise(tncnt_schisto, pvalue=0.9)
#x <- pickAboveTechNoise(ncnt_schisto, pvalue=0.001)
#x <- pickAboveTechNoise(ncnt_schisto, pvalue=0.00001)


#Check which genes got picked
x <- tncnt_schisto[notnoisy,]
x <- x[-grep("ERC",rownames(x)),]
pickedgenes <- sort(rownames(x))
pickedgenes
#nrow(x)



##########################################################################
################### Interesting gene sets ################################
##########################################################################


#tncnt is a bad way

###############################
## Define interesting subsets of genes
#rownames(tncnt_schisto)[grep("Ccr",rownames(tncnt_schisto))]
liCxcr <- rownames(tncnt_schisto)[grep("Cxcr",rownames(tncnt_schisto))]
liCcr <- rownames(tncnt_schisto)[grep("Ccr",rownames(tncnt_schisto))]

liSurf <- c("T-bet??","Gata3","Rora","Rorc","Runx3","AHR??","TGFbeta??","HLX??")
liCd <- rownames(tncnt_schisto)[grep("Cd",rownames(tncnt_schisto))]

liNeg <- c("Cd11b","Cd11c","Cd19","Cd8a","Ncam1","Cd33")
            #todo which Cd11??

liFoo <- c(liCxcr,"Rora",liCcr)


###############################
### Colors for different tissues and platxes
colorfortissue <- rep("black",ncol(cnt_schisto))
colorfortissue[grep("MLN",colnames(cnt_schisto))] <- "red"
colorfortissue[grep("SPL_",colnames(cnt_schisto))] <- "blue"
colorfortissue[grep("SPL2_",colnames(cnt_schisto))] <- "green"

colorforplate <- rep("black",ncol(cnt_schisto))
colorforplate[grep("SPL2_",colnames(cnt_schisto))] <- "green"






##########################################################################
################### Specific questions ###################################
##########################################################################


#how many cells express rora?
mean(log10(1+as.double(ncnt_schisto[toensid("Rora"),])) > gene_cutoff[["Rora"]])

mean(log10(1+as.double(ncnt_schisto[toensid("Rora"),
                                    cellcondition$isMLN])) > gene_cutoff[["Rora"]])

mean(log10(1+as.double(ncnt_schisto[toensid("Rora"),
                                    cellcondition$isSPL])) > gene_cutoff[["Rora"]])

#after log10(1+...)
gene_cutoff[["Rora"]] <- 0.5

ncnt_schisto[toensid()]
#cellcondition$


##Of the cytokine expressing cells, how many rora?
mean(as.double(ncnt_schisto[
  toensid("Rora"),
  cellcondition$isCD4 & cellcondition$highCytokine]) > gene_cutoff$Rora)   #71%







########################################################################################################
########################################################################################################
########################################################################################################
############################### Plots and stuff below  ---- no longer used
########################################################################################################
########################################################################################################
########################################################################################################


hist(log(tncnt_schisto[which(rownames(tncnt_schisto)=="Cd11c"),]+0.001),breaks=20)



#y <- x[rowMedians(x)>1,]
#y <- x[which(rownames(x) %in% liChemo),]
y <- tncnt_schisto[which(rownames(tncnt_schisto) %in% liFoo),]




#keepcells <- which(log(tncnt_schisto[which(rownames(tncnt_schisto)=="Sell"),])>3)
keepcells <- 1:ncol(tncnt_schisto)
heatmap.2(
  normalizegenes(log(rednovar(y)+0.1))[,keepcells],
  #log(y[,keepcells]+0.1),
  
  #  ColSideColors = colorforplate[keepcells],
  ColSideColors = colorfortissue[keepcells],
  
  key=TRUE,
  key.title="Expression",
  key.xlab="",
  key.ylab="",
  trace="none"
)




corgenered(y)





onegene <- tncnt_schisto[which(rownames(tncnt_schisto)=="Rora"),]
onegene2 <- tncnt_schisto[which(rownames(tncnt_schisto)=="Cxcr6"),]
cor(onegene, onegene2, method="spearman")

thecor <- cor(onegene, t(tncnt_schisto), method="spearman")
hist(thecor,breaks=30)
names(thecor) <- rownames(tncnt_schisto)
as.matrix(thecor[order(abs(thecor),decreasing = TRUE)][1:800])

# check subpopulations:
cell_cor =cor(tncnt_schisto, method = "spearman")

colorfortissue <- rep("black",ncol(cell_cor))
colorfortissue[grep("MLN",colnames(cell_cor))] <- "red"
colorfortissue[grep("SPL_",colnames(cell_cor))] <- "blue"
colorfortissue[grep("SPL2_",colnames(cell_cor))] <- "blue"
hmap=heatmap.2(cell_cor,trace="none",density="none",#scale="row",
               cexRow=0.5,cexCol=0.5,  RowSideColors = colorfortissue,
               col=redgreen, hclustfun=function(x) hclust(x,method="complete"))

color_fht <- rep("black",ncol(cell_cor))
color_fht[which(log(1+tncnt_schisto["Rora",]) > 2 )] = "pink"   
heatmap.2(cell_cor,trace="none",density="none",#scale="row",
          cexRow=0.5,cexCol=1.5,  ColSideColors = color_fht,
          col="redgreen", hclustfun=function(x) hclust(x,method="complete"))


colorforplate <- rep("black",ncol(cnt_schisto))
colorforplate[grep("SPL2_",colnames(cnt_schisto))] <- "green"

roraall=subset(All, All[,"RORA"] > -9)
rc = rep("grey",nrow(roraall))
rc[grep("G",rownames(roraall))]="orange"
hmap=heatmap.2(roraall,trace="none",density="none",#scale="row",
               cexRow=0.5,cexCol=0.5,  RowSideColors = rc,
               col=redgreen, hclustfun=function(x) hclust(x,method="complete"))

plot(tncnt_schisto["Rora",])
hist(log(1+tncnt_schisto["Rora",]),breaks=40)
hist(log(1+CD34cellsNor["ENSMUSG00000032238",]),breaks=40)
hist(log(1+tncnt_schisto["Foxp3",]),breaks=40)
plot(tncnt_schisto["Gata3",],tncnt_schisto["Il4",])

plot(log(1+tncnt_schisto["Foxp3",]),log(1+tncnt_schisto["Rora",]))
text(log(1+tncnt_schisto["Foxp3",]),log(1+tncnt_schisto["Rora",]),labels=colnames(tncnt_schisto),cex=0.5)

plot(log(1+tncnt_schisto["Foxp3",]),log(1+tncnt_schisto["Cxcr6",]))
text(log(1+tncnt_schisto["Foxp3",]),log(1+tncnt_schisto["Cxcr6",]),labels=colnames(tncnt_schisto),cex=0.5)

plot(log(1+tncnt_schisto["Foxp3",]),log(1+tncnt_schisto["Il2ra",]))
text(log(1+tncnt_schisto["Foxp3",]),log(1+tncnt_schisto["Il2ra",]),labels=colnames(tncnt_schisto),cex=0.5)


forpca <- tncnt_schisto
# apply(forpca,1,mean)
forpca <- tncnt_schisto[notnoisy,]
forpca <- forpca[-grep("ERC",rownames(forpca)),]

pcaRes = prcomp(t(normalizegenesSD(log(1+forpca))))
barplot(pcaRes$sdev)
plot(pcaRes$x[,1], pcaRes$x[,3], xlab="PC 1", ylab="PC 2", pch=19, cex=0.5)
pairs(pcaRes$x[,1:5])
colnames(pcaRes$rotation)

do.call(paste,as.list(names(sort(abs(pcaRes$rotation[,1]),decreasing=TRUE)[1:200])))
do.call(paste,as.list(names(sort(abs(pcaRes$rotation[,3]),decreasing=TRUE)[1:200])))
sort(abs(pcaRes$rotation[,2]),decreasing=TRUE)[1:20]


#######check how are my reads:

MT_genes = read.table("MT_genes.txt", sep="\t",head=T, row.names=1)
mit_genes <- cnt_schisto[rownames(MT_genes), ]
mouse_genes <- cnt_schisto[substr(rownames(cnt_schisto), 1, 5) == "ENSMU", ]
ercc_genes <- cnt_schisto[substr(rownames(cnt_schisto), 1, 4) == "ERCC", ]

all_reads  <- colSums(cnt_schisto)
mapped_reads <- colSums(mouse_genes)
ercc_reads <- colSums(ercc_genes)
mit_reads <- colSums(mit_genes)

per_exon_mapped_reads  <- mapped_reads/all_reads*100
per_ercc_mapped_reads <- ercc_reads/all_reads*100
per_unmapped_reads  <- cnt_schisto["__not_aligned",]/all_reads*100
per_mapped_not_exon_reads  <- cnt_schisto["__no_feature",]/all_reads*100
per_mito_mapped_reads  <- mit_reads/all_reads*100

crap <- cnt_schisto[substr(rownames(cnt_schisto), 1, 2) == "__", ]
crap_reads  <- colSums(crap)
percentage_crap_reads  <- crap_reads/all_reads*100

all_percentages = rbind(percentage_mito_mapped_reads, percentage_exon_mapped_reads,percentage_ercc_mapped_reads, percentage_unmapped_reads, percentage_mapped_not_exon_reads)
hist(percentage_mito_mapped_reads)
########remove cells 

#where mit reads and more then 10%
low_mit = names(mit_genes)[which(colSums(mit_genes)/all_reads*100<10)]

# where number of mapped reads are less then 100000 reads
high_mapped_reads = names(mouse_genes)[which(colSums(mouse_genes)>100000)]

# where per of unmapped is more then 30% 
low_per_unpammed = names(mit_genes)[which(cnt_schisto["__not_aligned",]/all_reads*100<30)]


super_cells = intersect(intersect(low_mit,high_mapped_reads),low_per_unpammed)
good_cells = cnt_schisto[,super_cells]

######### sort cell subtype

### pick up CD3/CD4 cells
CD3_CD4genes = c("ENSMUSG00000032093", "ENSMUSG00000023274", "ENSMUSG00000053977")

CD3_CD4t = ncnt_schisto[CD3_CD4genes,]
CD3_CD4t = t(CD3_CD4t)

#remove cells which CD3 and CD4 are not expressed
good_cells <- good_cells[-grep("_",rownames(good_cells)),] #remove rows with __
tgood_cells =t(good_cells)
CD34cells=subset(tgood_cells, tgood_cells[,"ENSMUSG00000032093"] > 10)
CD34cells = t(CD34cells)
CD34cells <- CD34cells[-grep("ERC",rownames(CD34cells)),]
CD34cellsNor= normalizegenesSD(log(1+CD34cells))

######## check which TF are expressed

key_genes = c("ENSMUSG00000015619", "ENSMUSG00000032238", "ENSMUSG00000039521", 
              "ENSMUSG00000001444", "ENSMUSG00000028150")
cell_key = log(1+CD34cells[key_genes,])

cell_key_counts = (cell_key>2)*1
sum_cell_key= apply(cell_key_counts, 1, sum)

#####    Pick genes which are highly exppresed in at least one cell

TPMfiltered <- CD34cellsNor[which(apply(CD34cellsNor, 1, max) >5 ), ]

#Spearman correlation matrix for cell pairs
cm <- cor(TPMfiltered, method="spearman")
color_fht <- rep("black",ncol(cm))
color_fht[which(log(1+TPMfiltered["ENSMUSG00000015619",]) > 5 )] = "pink" 
hmap=heatmap.2(cm,trace="none",density="none", ColSideColors = color_fht,
               cexRow=0.5,cexCol=0.5, col=redgreen)

# PCA
pcaRes <- prcomp(cm)
plot(pcaRes$x[,1], pcaRes$x[,2], xlab="PC 1", ylab="PC 2")
plot(pcaRes$x[,1], pcaRes$x[,3], xlab="PC 1", ylab="PC 2", pch=19, cex=0.5)
pairs(pcaRes$x[,1:5])

### Find cell type:
GATA3
plot(CD34cells["ENSMUSG00000015619",])
plot(log(1+CD34cells["ENSMUSG00000015619",]))
plot(log(1+CD34cellsNor["ENSMUSG00000015619",]))

RORA
plot(CD34cells["ENSMUSG00000032238",])
plot(log(1+CD34cells["ENSMUSG00000032238",]))
plot(log(1+CD34cellsNor["ENSMUSG00000032238",]))

FOXP3
plot(CD34cells["ENSMUSG00000039521",])
plot(log(1+CD34cells["ENSMUSG00000039521",]))
plot(log(1+CD34cellsNor["ENSMUSG00000039521",]))

Tbx21
plot(CD34cells["ENSMUSG00000001444",])
plot(log(1+CD34cells["ENSMUSG00000001444",]))
plot(log(1+CD34cellsNor["ENSMUSG00000001444",]))

RORC
plot(CD34cells["ENSMUSG00000028150",])
plot(log(1+CD34cells["ENSMUSG00000028150",]))
plot(log(1+CD34cellsNor["ENSMUSG00000028150",]))


IFNg
plot(CD34cells["ENSMUSG00000055170",])
plot(log(1+CD34cells["ENSMUSG00000055170",]))
plot(log(1+CD34cellsNor["ENSMUSG00000055170",]))


colnames(CD34cells)
cellcat <- rep("notype",ncol(CD34cells))
cellcat[log(CD34cellsNor["ENSMUSG00000015619",]+1)>2] <- "GATA3"
cellcat
CD34cells


#Venn diagram, assuming the matrix is normalized read counts
cutoff = 0.1
GATA3 <- (CD34cellsNor["ENSMUSG00000015619",]>cutoff)
RORA <- (CD34cellsNor["ENSMUSG00000032238",]>cutoff)
FOXP3 <- (CD34cellsNor["ENSMUSG00000039521",]>cutoff)
TBX21 <- (CD34cellsNor["ENSMUSG00000001444",]>cutoff)
RORC <- (CD34cellsNor["ENSMUSG00000028150",]>cutoff)

vd <- vennCounts(cbind(GATA3, RORA, FOXP3,TBX21, RORC))
vennDiagram(vd,counts.col = "red",names=c("GATA3", "RORA", "FOXP3","TBX21", "RORC"),cex=0.8)



