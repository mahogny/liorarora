source("common_functions.R")

if(FALSE){
  BiocManager::install("scater")  
  BiocManager::install("topGO")
}

library(mixtools)
library(DESeq2)
library(sqldf)
library(matrixStats)
library(statmod)
library(stringr)
library(gplots)
library(limma)
library(scater)


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




# cnt_file1 <- readrawcountmulti("schisto_c1/set1","schisto_c1/mappingSet1.csv")
# cnt_file2 <- readrawcountmulti("schisto_c1/set2","schisto_c1/mappingSet2.csv")
# cnt_schisto <- cbind(cnt_file1,cnt_file2)
# all(rownames(cnt_file1)==rownames(cnt_file2))  #Sanity check
# ncol(cnt_schisto)
# write.table(cnt_schisto, "schisto_c1/count.txt")

cnt_schisto <- read.table("schisto_c1/count.txt")


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
gene_count <- as.integer(colSums(dat[grep("ENSMUSG",rownames(dat)),]))
#Number of detected genes
detected_genes <- as.integer(colSums(dat>0))
#Proportion of mitochondrial reads
mt_counts <- as.integer(colSums(dat[mt_genes$ensembl_gene_id,]))
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
x<-findbadcells()
#as.double(gene_count[x]) #read counts of dropped genes
cellcondition$isgood[x]<-FALSE
dev.off()

### Remove cells with too little coverage as these will just disturb later statistics
#cellcondition$isgood[apply(cnt_schisto,2,sum)<1.5e6] <- FALSE




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



gene_cutoff<-list()

hist(log10(1+as.double(ncnt_schisto[toensid("Cd4"),])),breaks=20)
hist(log10(1+as.double(ncnt_schisto[toensid("Rora"),])),breaks=20, xlab="Log10(1+Rora)")
hist(log10(1+as.double(ncnt_schisto[toensid("Foxp3"),])),breaks=20)
hist(log10(1+as.double(ncnt_schisto[toensid("Tbx21"),])),breaks=20)

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
  apply(ncnt_schisto[toensid(str_to_title(heatmap_markers_cyto)) %in% rownames(ncnt_schisto),]> 0.5,2,any)
# log10(1+as.double(ncnt_schisto[toensid("Il4"),])) > 0.5 |
#   log10(1+as.double(ncnt_schisto[toensid("Il10"),])) > 0.5 |
#   log10(1+as.double(ncnt_schisto[toensid("Il13"),])) > 0.5

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
colorfortissue[grep("MLN",  colnames(cnt_schisto))] <- "red"
colorfortissue[grep("SPL_", colnames(cnt_schisto))] <- "blue"
colorfortissue[grep("SPL2_",colnames(cnt_schisto))] <- "green"

colorforplate <- rep("black",ncol(cnt_schisto))
colorforplate[grep("SPL2_",colnames(cnt_schisto))] <- "green"




##########################################################################
################### RNAseq QC and normalization ##########################
##########################################################################

### diagram, expression of a gene over time and tissue?

#cellcondition$

#there is MLN and SPL!


##########################################################################
################### Gene exp heatmap for paper ###########################
##########################################################################

heatmap_markers_cyto <- str_to_upper(c())
heatmap_markers <- str_to_upper(c(
  "Rora", "Gata3","Tbx21", "Foxp3", "Pparg",
  "Cd44", "Sell", "Icos",
  "Cxcr6", "Ccr2", "Cxcr5", "Il1rl1" ##bonus
))   

col_spleen <- "red"  #hacks, not what in paper
col_mesLN <- "blue"


celltypecol <- rep("black", ncol(ncnt_schisto))
celltypecol[cellcondition$isSPL]  <- col_spleen
celltypecol[cellcondition$isMLN]  <- col_mesLN

s <- t(ncnt_schisto[rownames(ncnt_schisto) %in% toensid(str_to_title(c(heatmap_markers,heatmap_markers_cyto))),])
so <- 0 + (s>0.5)
usecell <-  cellcondition$isCD4 #!cellcondition$highSell &

so <- data.frame(so[usecell,toensid(str_to_title(heatmap_markers))])
colnames(so) <- str_to_title(togenesym(colnames(so)))

celltypecol_sub <- celltypecol[usecell]

so_order <- order(so$Rora, so$Gata3, so$Foxp3, so$Tbx21, so$Pparg, so$Cd44, so$Sell,
                  so$Icos, so$Cxcr6, so$Ccr2, so$Cxcr5, so$Il1rl1)#, so$Cytokines)

########## PDF generation
pdf("out.schisto/heatmap.pdf",h=4)
heatmap.2(
  t(so[so_order,]), 
  col=c("#CCCCCC","red"),
  rowsep = 1:ncol(s),
  trace = "none",
  scale="none",
  Rowv = FALSE,
  Colv = FALSE,
  cexRow = 1.5,
  ColSideColors = celltypecol_sub[so_order], 
  dendrogram = "none",
  labCol = FALSE)
dev.off()


###### Show a representative gene (reviewer question)



##########################################################################
################### Specific questions ###################################
##########################################################################

ensidGata3 <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Gata3"]
ensidRora <- toensid("Rora")
ensidFoxp3 <- toensid("Foxp3")


mean(ncnt_schisto[ensidRora,] > 0.5)  #59%


# fisher.test(
#   (ncnt_schisto[ensidFoxp3,] > 0.5)[cellcondition$isgood],                  
#   (ncnt_schisto[ensidRora,] > gene_cutoff[["Foxp3"]])[cellcondition$isgood])
# cor.test(ncnt_schisto[ensidFoxp3,cellcondition$isgood], ncnt_schisto[ensidRora,cellcondition$isgood])

################## Overlap Rora & Foxp3?
overl <- data.frame(
  rora =ncnt_schisto[ensidRora,] > 0.5,
  foxp3=ncnt_schisto[ensidFoxp3,] > 0.5,
  gata3=ncnt_schisto[ensidGata3,] > 0.5)
overl <- overl[cellcondition$isgood,]

overl_foxp3 <- overl[!overl$gata3, c("rora","foxp3")]
overl_foxp3

overl_foxp3_red <- overl[overl$foxp3 | overl$rora, c("rora","foxp3")]
fisher.test(overl_foxp3_red$rora, overl_foxp3_red$foxp3)   # 0.0001918 
sum(overl$foxp3 & overl$rora)  #16
sum(overl$rora) # 49

# overl_gata3 <- overl[!overl$foxp3, c("rora","gata3")]
# overl_gata3
# 
# fisher.test(overl$rora[!overl$foxp3], overl$gata3[!overl$foxp3])

# overl <- data.frame(
#   rora =ncnt_schisto[ensidRora,] > 0.5,
#   gata3=ncnt_schisto[ensidGata3,] > 0.5)
# 
# vd <- vennCounts(overl[!overl$foxp3,c(1,3)])
# vennDiagram(vd)

#,counts.col = "red",names=c("GATA3", "RORA", "FOXP3","TBX21", "RORC"),cex=0.8)


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

