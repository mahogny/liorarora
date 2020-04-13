library(stringr)
library(rARPACK)
# library(rgl)
#library(plot3D)
library(scatterplot3d)
library(roots)
library(reshape2)
library(gplots)
library(mixtools)



# set the directory to your current working directory (usually where your data files are)
#setwd("/Users/lv5/Documents/Lab/Sarah/Projects/Single cell gene expression/7. Nippo exp/Analysis/Liora")



############################################################################################
##################### read biomark data and clean up #######################################
############################################################################################



# read in the data
tbl_1 <- read.csv (file="nipo_biomark/exp_1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
tbl_2 <- read.csv (file="nipo_biomark/exp_2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# update the values
# if on any row column 4 value is 'Fail', then column 3 value will be set as 999
tbl_1[tbl_1[, 4]=='Fail', 3] <- 999
tbl_2[tbl_2[, 4]=='Fail', 3] <- 999

# remove the call column
tbl_1 <- tbl_1[,1:3]
tbl_2 <- tbl_2[,1:3]


#reshaping the data file so that rows are samples and columns are genes
tw1 <- acast (data=tbl_1, Cell ~ Gene, value.var="Ct")
tw2 <- acast (data=tbl_2, Cell ~ Gene, value.var="Ct")
identical(colnames(tw1), colnames(tw2))
tw0 <- rbind(tw1, tw2)

#change 999 to 40
tw0[tw0[,]==999] <- 40

#Remove cells in which Ubc is not expressed 
tw1 <- subset(tw0, tw0[,"UBC"] < 26)
tw1 <- subset(tw0,tw0[,"UBC"] < 26 | tw0[,"PTPRC"] < 26)

#Remove cells in which negative markers are expressed
tw2 <- subset(tw1, tw1[,"ITGAX"] == 40)
tw2 <- subset(tw2, tw2[,"ITGAM"] == 40)
tw2 <- subset(tw2, tw2[,"CD19"] == 40)
tw2 <- subset(tw2, tw2[,"CD8A"] == 40)
tw2 <- subset(tw2, tw2[,"NCAM1"] == 40)
tw2 <- subset(tw2, tw2[,"CD33"] == 40)

#Remove cells in which negative markers are expressed
tw2<-tw2[rowSums(tw2)!=40*ncol(tw2),colSums(tw2)!=40*nrow(tw2)] 


############### Normalization
normdata <- matrix(0, nrow(tw2), ncol(tw2))
rownames(normdata)<-rownames(tw2)
colnames(normdata)<-colnames(tw2)

listrefgene <- c('ATP5A1', 'HPRT1', 'UBC')
for(i in 1:nrow(tw2)){
  #reference genes
  hk_ct <- c(tw2[i,'ATP5A1'], tw2[i,'HPRT1'], tw2[i,'UBC'])
  
  #check if any of the reference points are crap. if so, remove
  toremove <- which(hk_ct >= 27)
  if (length(toremove) > 0 & length(toremove) < 3){
    hk_ct <- hk_ct[-toremove]
  }
  
  #normalize this cell
  for(j in 1:ncol(tw2)){
    normdata[i,j] <- 2^(mean(hk_ct) - tw2[i,j])
  }
}

############### Rename cells and sort by name
oldnames <- strsplit(rownames(normdata), "_")
newnames <- rep("",length(oldnames))
for (i in 1:618){
  newnames[i] <- sprintf("%s_%s_%s_%s_%s", oldnames[[i]][5], oldnames[[i]][1],oldnames[[i]][3],oldnames[[i]][2],oldnames[[i]][4])
}
rownames(normdata) <- newnames
#normdata <- normdata[order(rownames(normdata)),]
tblsorted <- normdata[order(rownames(normdata)),]


############################################################################################
######################### optimal cutoffs? #################################################
############################################################################################

#hist(log10(as.double(CD3_4proP[,"RORA"])))
#cutoff_rora <- -3



seed(0)
gene_cutoff<-list()
for(c in colnames(tblsorted)){
  #mean(normalmixEM(log10(as.double(tblsorted[,"RORA"])))$mu)
  gene_cutoff[[c]] <- 10^( 1/2 * mean(normalmixEM(log10(as.double(tblsorted[,c])))$mu)  )
#  v <- log10(as.double(tblsorted[,c]))
#  gene_cutoff[[c]] <- 10^(min(v)*0.1 + max(v)*0.9)
}
gene_cutoff$RORA
gene_cutoff$CD4
gene_cutoff$CD3E

### loose cutoffs 
###CD4 <- subset(CD3_4proP, CD3_4proP[,"CD4"] > -9 | CD3_4proP[,"CD3E"] > -9)

hist(log10(tblsorted[,"RORA"]),xlab = "Log10(Rora)")

############################################################################################
######################### get subsets ######################################################
############################################################################################


###### Remove cells in which CD3 & CD4 marker is not expressed = Th cells

cellcondition <- data.frame(isgood=rep(TRUE,nrow(tblsorted)))
cellcondition$isTcell <- FALSE
cellcondition$isTcell[c(grep("ac_",rownames(tblsorted)), grep("nv_",rownames(tblsorted)))] <- TRUE

###### Assign tissue
cellcondition$isGut <- FALSE
cellcondition$isGut[grep("g_",rownames(tblsorted))] <- TRUE
cellcondition$isLung <- FALSE
cellcondition$isLung[grep("l_",rownames(tblsorted))] <- TRUE
cellcondition$isME <- FALSE
cellcondition$isME[grep("me",rownames(tblsorted))] <- TRUE
cellcondition$isMED <- FALSE
cellcondition$isMED[grep("med",rownames(tblsorted))] <- TRUE
cellcondition$isMES <- FALSE
cellcondition$isMES[grep("mes",rownames(tblsorted))] <- TRUE

cellcondition$tissue <- ""
cellcondition$tissue[cellcondition$isMES]  <- "mesLN"
cellcondition$tissue[cellcondition$isMED]  <- "medLN"
cellcondition$tissue[cellcondition$isLung] <- "lung"
cellcondition$tissue[cellcondition$isGut]  <- "gut"

cellcondition$highRora <- FALSE
cellcondition$highRora <- as.double(tblsorted[,"RORA"])> gene_cutoff$RORA
cellcondition$highIl4 <- FALSE
cellcondition$highIl4 <- as.double(tblsorted[,"IL4"])> gene_cutoff$IL4
cellcondition$highIl13 <- FALSE
cellcondition$highIl13 <- as.double(tblsorted[,"IL13"])> gene_cutoff$IL13

cellcondition$highSell <- FALSE
cellcondition$highSell <- as.double(tblsorted[,"SELL"])> 10^3

#### Assign day
cellcondition$day <- -1
cellcondition$day[grep("D0_",rownames(tblsorted))] <- 0
cellcondition$day[grep("D03",rownames(tblsorted))] <- 3
cellcondition$day[grep("D05",rownames(tblsorted))] <- 5
cellcondition$day[grep("D07",rownames(tblsorted))] <- 7
cellcondition$day[grep("D10",rownames(tblsorted))] <- 10

#cellcondition$isCD4 <- cellcondition$isTcell & (tblsorted[,"CD4"] > -9 | tblsorted[,"CD3E"] > -9)
cellcondition$isCD4 <- cellcondition$isTcell & 
  (tblsorted[,"CD4"] > gene_cutoff["CD4"] | tblsorted[,"CD3E"] > gene_cutoff["CD3E"])

cellcondition$highCytokine <- cellcondition$highIl4 | cellcondition$highIl13

################ maybe don't use this

#(naive and activated)
# By protein
# CD3_4proP <- tblsorted[c(grep("ac_",rownames(tblsorted)), grep("nv_",rownames(tblsorted))),]  #nv?
# # CD4 cell would be cell either expressing CD3 or CD4 (excludes CD8 before)
# CD4 <- subset(CD3_4proP, CD3_4proP[,"CD4"] > -9 | CD3_4proP[,"CD3E"] > -9)
# #write.csv(CD4,"CD4.csv")



############################################################################################
##################### how much rora exp in different populations? ##########################
############################################################################################

#How many cells are rora in different tissues?

mean(as.double(tblsorted[
  cellcondition$isGut | cellcondition$isLung,
  "RORA"])> gene_cutoff$RORA)   #37%

mean(as.double(tblsorted[
  cellcondition$isME,
  "RORA"])> gene_cutoff$RORA)   #20%

##Of the cytokine expressing cells, how many rora?
mean(as.double(tblsorted[
  cellcondition$isCD4 & cellcondition$highCytokine,
  "RORA"])> gene_cutoff$RORA)   #83%

#cellcondition$highCytokine


############################################################################################
######################### Gene exp heatmap for paper #######################################
############################################################################################


heatmap_markers_cyto <- str_to_upper(c("Ifng", "Il4", "Il13", "Il5", "Il6", "Il10", "Il17a"))
heatmap_markers <- str_to_upper(c(
  "Rora", "Gata3","Tbx21", "Foxp3", "Pparg",
  "Cd44", "Sell", "Icos",
  "Cxcr6", "Ccr2", "Cxcr5", "Il1rl1" ##bonus
  ))   

col_mesLN  <- "#0000FF"
col_medLN  <- "#5555FF"
col_lung   <- "#BB0000"
col_gut    <- "#FF5555"
col_spleen <- "yellow"  #additional
#col_medLN  <- "#8888FF"
#col_gut    <- "yellow"

celltypecol <- rep("black", ncol(tblsorted))
celltypecol[cellcondition$isGut]  <- col_gut
celltypecol[cellcondition$isLung] <- col_lung
celltypecol[cellcondition$isMES]  <- col_mesLN
celltypecol[cellcondition$isMED]  <- col_medLN

s <- tblsorted
for(i in 1:ncol(s)){
  s[,i] <- 0+(s[,i] >  1e-4)  #gene_cutoff[[colnames(s)[i]]])
}
usecell <-  cellcondition$isCD4 & apply(s,1,sum)<50  #!cellcondition$highSell &
scyto <- (apply(s[,heatmap_markers_cyto],1,sum)>0)+0

s <- data.frame(s[,heatmap_markers], Cytokines=scyto)
#s <- s[usecell,c("Cytokine",heatmap_markers)]
colnames(s) <- str_to_title(colnames(s))  

so <- s[usecell,]
so_order <- order(so$Rora, so$Gata3, so$Foxp3, so$Tbx21, so$Pparg, so$Cd44, so$Sell,
                  so$Icos, so$Cxcr6, so$Ccr2, so$Cxcr5, so$Il1rl1, so$Cytokines)
so <- so[so_order,]
# so <- so[order(celltypecol[usecell],so$Rora, so$Gata3, so$Foxp3, so$Tbx21, so$Pparg, so$Cd44, so$Sell,
#                so$Icos, so$Cxcr6, so$Ccr2, so$Cxcr5, so$Il1rl1, so$Cytokines),]


pdf("out.nipobiomark/heatmap.pdf",h=4)  ####### need to reimport!!
heatmap.2(
  t(so), 
  col=c("#CCCCCC","red"),
  rowsep = 1:ncol(s),
  trace = "none",
  scale="none",
  Rowv = FALSE,
  Colv = FALSE,
  cexRow = 1.5,
  ColSideColors = celltypecol[usecell][so_order],  #modified . celltypecol is wrong!
  dendrogram = "none",
  labCol = FALSE)
dev.off()


getmatrix_daytissue <- function(genename="Rora", usecell=rep(TRUE, nrow(s)),
                                checktissue=c("mesLN","medLN","gut","lung"),
                                checkday=c(3,5,7), getcount=FALSE){
  out <- matrix(nrow=length(checktissue), ncol=length(checkday))
  for(curtis in 1:length(checktissue)){
    for(curday in 1:length(checkday)){
      v <- cellcondition$tissue==checktissue[curtis] & cellcondition$day==checkday[curday] & usecell
      #v <- v[usecell]
      if(getcount){
        #num cells to check
        out[curtis,curday] <- sum(v) 
      } else {
        #fraction of cells expressing
        out[curtis,curday] <- mean(s[v, genename])  
      }
#      print(sprintf("%s   %s   %s", checktissue[curtis], checkday[curday],  out[curtis,curday]))
    }
  }
  colnames(out) <- sprintf("d%s",checkday)
  rownames(out) <- c("Lung","Gut","MedLN","MesLN")
  out
}

plotonegenescatter_daytissue <- function(genename="Rora", usecell=rep(TRUE, nrow(s)), fcex=2.5,
                                         checktissue=c("mesLN","medLN","gut","lung"),
                                         checkday=c(3,5,7)){

  #Restrict to tissues and times  
  usecell <- usecell & 
    cellcondition$tissue %in% checktissue & 
    cellcondition$day    %in% checkday
  
  #s[cellcondition$tissue=="medLN" & cellcondition$day==5,"Sell"]
  
  #Calculate percentages
  out <- getmatrix_daytissue(genename, usecell=usecell)
  cnt <- getmatrix_daytissue(genename, usecell=usecell, getcount = TRUE)
  
  print(out)
  print(cnt)
  
  #Calculate XY positions
  set.seed(0)
  cellx <- rep(-1, nrow(s))
  celly <- rep(-1, nrow(s))
  for(curtis in 1:length(checktissue)){
    for(curday in 1:length(checkday)){
      v <- cellcondition$tissue==checktissue[curtis] & cellcondition$day==checkday[curday]
      cellx[which(v)] <- curday + runif(n=sum(v),min = 0.1, max=0.9) - 1
      celly[which(v)] <- curtis + runif(n=sum(v),min = 0.1, max=0.9) - 1
    }
  }

  #Make the grid
  plot(x=c(0,0),y=c(0,0),type="l",xlim=c(-2,6), ylim=c(-2,6), axes = FALSE)
  for(i in 0:nrow(out)){
    lines(x=c(0,ncol(out)),y=c(i,i),type="l")
  }
  for(i in 0:ncol(out)){
    lines(x=c(i,i),y=c(0,nrow(out)),type="l")
  }

  #Gene name label
  text(x=-1,y=nrow(out)+1, labels=genename, cex=fcex*1.5)
  
  #Axis labels
  text(x=-0.8,y=(1:nrow(out))-0.5,labels=rownames(out),cex=fcex)
  text(x=(1:ncol(out))-0.5,y=nrow(out)+0.3,labels=colnames(out),cex=fcex)
  
  #Percentages to the right
  text(x=ncol(out)+0.5, y=1:nrow(out)-0.5,cex=fcex, 
       labels=sprintf("%s%%",round(100*apply(out,1,mean))))
  #print(apply(out,1,mean))
  
  #Percentages beneath
  text(y=-0.3, x=1:ncol(out)-0.5,cex=fcex, 
       labels=sprintf("%s%%",round(100*apply(out,2,mean))))
  #print(apply(out,2,mean))
  
  #Assign colors to points
  pp <- usecell  #cellx>=0 & 
  thecol <- s[,genename]
 # print(thecol)
  thecol[thecol==0] <- "gray"
  thecol[thecol==1] <- "red"
#  print(thecol)

  points(x=cellx[pp], y=celly[pp], pch=19, col=thecol[pp])
}


#getmatrix_daytissue("Sell")
plotonegenescatter_daytissue("Sell")

#Make multiple XY plots. Reuse XY for consistency
plotmultiscattergenes <- function(genestoplot=c("Rora","Sell","Gata3","Tbx21","Foxp3","Cytokines")){
  usecell <- rep(TRUE, nrow(s))
  usecell <- s$Rora==1
  #usecell <- s$Sell==1
  pdf("out.nipobiomark/gene_timespace.pdf",w=length(genestoplot)*7)
  par(mfrow=c(1,length(genestoplot)))
  for(g in genestoplot){
    print(g)
    plotonegenescatter_daytissue(g, usecell = usecell)
  }
  dev.off()
}
plotmultiscattergenes()



############################################################################################
######################### overlap between rora and other genes #############################
############################################################################################

b <- tblsorted[cellcondition$isCD4,]
b <- data.frame(
  Rora =b[,"RORA" ] > gene_cutoff$RORA,
  Sell=b[,"SELL"] > gene_cutoff$SELL
)
table(b$Rora[!b$Sell])  ##61/(94+61)
vc <- vennCounts(b)



b <- tblsorted[!cellcondition$highSell & cellcondition$isCD4,]
b <- data.frame(
  Gata3=b[,"GATA3"] > gene_cutoff$GATA3,
  Rora =b[,"RORA" ] > gene_cutoff$RORA,
  Foxp3=b[,"FOXP3"] > gene_cutoff$FOXP3,
)
vc <- vennCounts(b)
vennDiagram(vc)#,cex=c(1.5,1.5,1.5))


### Calculate p-value of rora & gata3 overlap being by random chance. Fisher test
fisher.test(b$Gata3, b$Rora)  #0.0003884
fisher.test(b$Foxp3, b$Rora)  #0.0001945


#only consider active... sell low!

hist(log10(1+as.double(tblsorted[
  cellcondition$isCD4 & cellcondition$highCytokine,
  "SELL"])))

mean(as.double(tblsorted[
  cellcondition$isCD4 & cellcondition$highCytokine,
  "RORA"])> gene_cutoff$RORA)   #83%



cor.rora  <- sort(decreasing = TRUE,cor(tblsorted,method = "spearman")["RORA",])
cor.gata3 <- sort(decreasing = TRUE,cor(tblsorted,method = "spearman")["GATA3",])
cor.foxp3 <- sort(decreasing = TRUE,cor(tblsorted,method = "spearman")["FOXP3",])

#b <- 
#apply(tblsorted,1,sum)

#cor.rora

b <- tblsorted[
  cellcondition$isCD4 & cellcondition$isTcell & cellcondition$isgood,
  union(union(names(cor.gata3[1:5]),names(cor.foxp3[1:5])),names(cor.rora[1:10]))]
for(i in 1:ncol(b)){
  b[,i] <- b[,i]>gene_cutoff[[colnames(b)[i]]] #  /max(b[,i])
}
#b <- log10(1+b)
#b <- b[apply(b,1,sum)>3,]  #I think this is a bit cheating

sort(apply(b,1,sum))


plot.new()
heatmap.2(
  b, 
  labRow = FALSE,
  trace="none",
  dendrogram="none") 
#heatmap.2(b)


############################################################################################
######################### diffusion map ####################################################
############################################################################################



diffMapFlor <- function(data, ndims = 4, nn=0.2, sigma=12, removeFirst = TRUE) {
  nn <- ceiling(nrow(data) * nn)      # Number of nearest neighbours to include
  KT <- sigma^2                       # Diffusion scale parameter
  cat("Calculating distance matrix.\n")
  d2 <- as.matrix(dist(data))^2       # distance matrix calculation
  R <- apply(d2, 2, function(x) sort(x)[nn])
  cat("Calculate nearest neighbours.\n")
  R <- matrix(rep(R,ncol(d2)), ncol = ncol(d2)) # Find distance for nn closest neighbour
  W <- exp(-d2 / (2*KT))              # Apply Gaussian kernel
  W <- (d2<R) * W                     # Only keep nn nearest neighbours
  D <- colSums(W)
  cat("Calculating and applying local density correction.\n")
  q <- D %*% t(D)                     # Calculate local density
  diag(W) <- 0
  H <- W / q                          # Correct for local density
  colS <- colSums(H)
  eS <- (colS == 0)                   # Excluded cells with no transitions
  Hp <- H[!eS,!eS] / colS[!eS]        # Normalise matrix
  cat("Calculating eigen vectors. Please wait...")
  n <- nrow(d2)
  cat(ndims+1)
  decomp <- eigs(Hp, which = "LR", ndims + 1)
  ## cat("Done.\n Calculating old way...")
  ## E <- eigen(Hp)                      # Eigen decomposition
  cat("Done.\n")
  if (removeFirst)
    startDim <- 2
  else 
    startDim <- 1
  chooseDims <- function(E) {
    eigOrd <- order(Re(E$values), decreasing = TRUE)
    E$values <- Re(E$values[eigOrd][startDim:(ndims+1)])
    E$vectors <- Re(E$vectors[,eigOrd][,startDim:(ndims + 1)]) # Remove first eigen vector
    rownames(E$vectors) <- rownames(data)[!eS]
    colnames(E$vectors) <- 1:ncol(E$vectors)
    return(E)
  }
  ## E <- chooseDims(E)
  decomp <- chooseDims(decomp)
  if (length(which(eS)) != 0) 
    cat(paste("\nWarning:\nCells \"", paste(names(which(eS)), collapse = ", "), "\" have no transitions and have been removed from the analysis\n", sep=""))
  return(decomp)
}


data <- tblsorted
data <- tblsorted[cellcondition$isCD4,]  #not converging
#data <- CD4                         #TODO is this the right set to use?
d <- diffMapFlor(data,ndims=3)
#d <- diffMapFlor(data,ndims=2)


#rc <- rep("brown",nrow(CD4set))
# data=CD4set #reg#early.activated#nuocyte#liora.data.all
# gene="T helper cells" #"UBC"#"RORA"#"UBC"
#rc[grep("me",rownames(CD4set))] <- "white"

# data=CD4set #reg#early.activated#nuocyte#liora.data.all
# gene="RORA"
# rc[as.double(data[,"RORA"]) > -5] <- "gray"

#plot3d(x=d$vectors[,1], y=d$vectors[,2], z=d$vectors[,3], xlab="PC 1", ylab="PC 2",zlab="PC 3",col=rc, main=gene, size=0.5, type="s")
#rgl.postscript("out/diffmap_1.eps")
#rgl.postscript("out/diffmap_1.pdf","pdf")

dodiffplot3 <- function(fname, main, rc){
  #pdf(fname)
  scatterplot3d(
    x=-d$vectors[,3], 
    y=d$vectors[,1], 
    z=-d$vectors[,2], 
    axis = TRUE,
    xlab="", ylab="",zlab="",
    #    xlab="PC 1", ylab="PC 2",zlab="PC 3",
    cex.symbols = 3,
    cex.lab = 1,
    cex.axis = 1,
    main=main, pch=3, color = rc)   #pch3
  #dev.off()
}

dodiffplot <- function(fname, main, rc){
  #pdf(fname)
  plot(
    x=-d$vectors[,1], 
    y=d$vectors[,2], 
    axis = TRUE,
    xlab="", ylab="",zlab="",
    #    xlab="PC 1", ylab="PC 2",zlab="PC 3",
    cex.symbols = 3,
    cex.lab = 1,
    cex.axis = 1,
    main=main, pch=3, color = rc)   #pch3
  #dev.off()
}


dodiffplotgene <- function(gene, col="brown", cutoff=-5){
  rc <- rep(col,nrow(data))
  rc[as.double(data[,toupper(gene)]) < cutoff] <- "gray"    
  dodiffplot(sprintf("out/diffmap_%s.pdf",gene), gene, rc)
}


pdf(file = "out.biomark/diffmap.pdf",w=5*7)
## set up the new plotting device (pdf)
par(mfrow = c(1,5))

dodiffplotgene("Rora",  "blue")
dodiffplotgene("Foxp3", "green")
dodiffplotgene("Gata3", "red")
dodiffplotgene("Sell",  "purple") 

rc <- rep("black",nrow(data))             # the rest: l and g (lung and gut)
rc[grep("me",rownames(data))] <- "gray"   # both type: mes and med
dodiffplot("out.biomark/diffmap_subset_me.pdf", "Lung+gut vs LN", rc)

dev.off()



############################################################################################
######################### another approach at clustering ###################################
############################################################################################


# fig 1b
# 
# cells: activated. 
# 
# heatmap: rora+ only, foxp3+ only   (only 2 rows in new)  , neither
# 
# cytokines: il4, il10, il13, il17, ifng, 
# 
# 
# which genes are co-exp rora & foxp3 together?
# 
# 
# which genes are DE?


############################################################################################
######################### redoing the clusterings ##########################################
############################################################################################

#Grab the RORA+ population
#pop_rora <- CD3_4proP
#data <- tblsorted[cellcondition$isCD4,]
pop_rora <- tblsorted[
  cellcondition$isCD4 & 
    tblsorted[,"RORA"]>gene_cutoff$RORA & 
    tblsorted[,"SELL"]<gene_cutoff$SELL,]

pop_foxp3 <- tblsorted[
  cellcondition$isCD4 & 
    tblsorted[,"FOXP3"]>gene_cutoff$RORA & 
    tblsorted[,"SELL"]<gene_cutoff$SELL,]





#normalize each gene
#pop_rora <- log10(pop_rora)
for(i in 1:ncol(pop_rora)){
  pop_rora[,i] <- pop_rora[,i]/max(pop_rora[,i])
}
#pop_rora

dim(pop_rora)
dim(CD3_4proP)

pop_rora[,c("GATA3","FOXP3","RORA","")]
#cytokines? wtf??

heatmap(log10(pop_rora))

ncol(CD3_4proP)

nrow(CD3_4proP)



############################################################################################
######################### more weird shit ################################################
############################################################################################



### Clamp lower value
CD4set <- CD4
CD4set[CD4set< -9] <- -9
write.csv(CD4set,"CD4set.csv")
hmap <- heatmap.2(CD4, col="redgreen", trace="none", cexRow=0.5, cexCol=0.5)
#plot(hmap$rowDendrogram)

pop1 = labels(cut(hmap$rowDendrogram, 150)$lower[[1]])
pop2 = labels(cut(hmap$rowDendrogram, 150)$lower[[2]])

CD4 = CD4[pop1,]

write.csv(CD4,"CD4.csv") ################# why twice ????????????????????????????????????????
CD4set = CD4

CD4set[CD4set< -9] <- -9
write.csv(CD4set,"CD4set.csv")


# heatmap for all cells
rc = rep("black",nrow(CD4set))
rc[grep("ac",rownames(CD4set))]="grey"
hmap=heatmap.2(CD4set,trace="none",density="none",
               cexRow=0.5,cexCol=0.5,  RowSideColors = rc,
               col=redgreen, hclustfun=function(x) hclust(x,method="complete"))

#Find the diffrences in cell sub-populations activated and naives
plot(hmap$rowDendrogram)

pop1 = labels(cut(hmap$rowDendrogram, 53.5)$lower[[1]])
pop2 = labels(cut(hmap$rowDendrogram, 53.5)$lower[[2]])
pop3 = labels(cut(hmap$rowDendrogram, 53.5)$lower[[3]])
pop2= c(pop2,pop3)
CD4p1 = CD4[pop1,]
CD4p2 = CD4[pop2,]

#ac/na diff between naive and activated:
CD4p1_2 <- sapply(seq(ncol(CD4p1)), function(x) f(CD4p1[,x], CD4p2[,x]))
CD4p1_2n = colnames(CD4p1)[CD4p1_2[3,] < 0.05]

CD4p1dif=CD4p1[,CD4p1_2n]
CD4p2dif=CD4p2[,CD4p1_2n]

CD4p1mean=apply(CD4p1dif, 2, mean)
CD4p2mean=apply(CD4p2dif, 2, mean)

diff= rbind(CD4p1mean, CD4p2mean)
diff=t(diff)
write.csv(CD4set, "CD4set1.csv")
CD4set2 = read.csv("CD4set2.csv", header=T, row.names=1, sep=",")
CD4set2= as.matrix(CD4set2)


hmap=heatmap.2(CD4set2,trace="none",density="none", Colv = FALSE,
               cexRow=0.5,cexCol=0.5,
               col=bluered)




# ########## find diff genes between activated...
rc = rep("black",nrow(CD4set))
rc[grep("ac",rownames(CD4set))]="grey"
hmap=heatmap.2(CD4set,trace="none",density="none",
               cexRow=0.5,cexCol=0.5,  RowSideColors = rc,
               col=redgreen, hclustfun=function(x) hclust(x,method="complete"))

#Find the diffrences in cell sub-populations activated and naives
#plot(hmap$rowDendrogram)

pop1 = labels(cut(hmap$rowDendrogram, 53.5)$lower[[1]])
pop2 = labels(cut(hmap$rowDendrogram, 53.5)$lower[[2]])
pop3 = labels(cut(hmap$rowDendrogram, 53.5)$lower[[3]])

CD4ap1 = CD4[pop1,]
CD4ap2 = CD4[pop2,]
CD4ap3 = CD4[pop3,]

color_fht = rep("white", 277)
color_fht[which(CD4ap1[,"FOXP3"] > -9 )] = "green"   

hmap=heatmap.2(CD4ap1,trace="none",density="none",
               cexRow=0.5,cexCol=0.5,  RowSideColors = color_fht,
               col=redgreen)


############################################################################################
######################### other stuff below ################################################
############################################################################################

#liora.data.all <- read.csv("Liora_data_byday.csv", header=T, row.names=1)
#liora.data.all[liora.data.all < -9] <- -9
#non.nuocyte=liora.data.all[readLines("non_nuocyte.txt"), ]
log.data <- CD4set
stage = sapply(rownames(log.data), function (r) { if (grepl("_ac", r)) { "ac" } else { "nv" } })
colour = sapply(rownames(log.data), function (r) { if (grepl("_ac", r)) { "red" } else { "black" } })
cellData = data.frame(cellID=rownames(log.data), embryoStage=stage, color=colour, stringsAsFactors = FALSE)
genes <- data.frame(colnames(log.data), stringsAsFactors = FALSE)
hscd <- new("SCD", experimentType = "qPCR", assayData=t(log.data), genoData=genes, phenoData=cellData)
hscd <- runPCA(hscd)
hscd <- diffuse(hscd, ndims=3)
temp <- roots:::selectCells(hscd, reduceMethod = "diffMap", dims = 3, useDims = c(1,2,3), doTree = FALSE, axLabs = "Diff Mode", selectedPath = "ESCpath", doPath=FALSE)
acp3=names(temp)

acp3pop=CD4[acp3,]
non.acp3pop <- CD4[setdiff(rownames(CD4), acp3), ]




##### functions ###############
#for :   ac/ na diff between naive and activated:
f <- function(x,y){
  test <- t.test(x,y, paired=F)
  data.frame(stat = test$statistic,
             df   = test$parameter,
             pval = test$p.value,
             conl = test$conf.int[1],
             conh = test$conf.int[2])
}




#ac/na:
p3 <- sapply(seq(ncol(acp3pop)), function(x) f(acp3pop[,x], non.acp3pop[,x]))
p3n = colnames(acp3pop)[p3[3,] < 0.05]

p3_diff=acp3pop[,p3n]
nonp3diff=non.acp3pop[,p3n]

p3mean=apply(p3_diff, 2, mean)
nonp3mean=apply(nonp3diff, 2, mean)

diff= rbind(p3mean, nonp3mean)
diff=t(diff)

hmap=heatmap.2(CD4set, col="redgreen", trace="none", cexRow=0.5, cexCol=0.5)

color_fht = rep("white", 420)
color_fht[which(CD4set[,"FOXP3"] > -9 )] = "green"   

rc = rep("white",nrow(CD4set))

rc[grep("D03_",rownames(CD4set))]="grey"
rc[grep("D05_",rownames(CD4set))]="yellow"
rc[grep("D07_",rownames(CD4set))]="red"

rc = rep("white",nrow(CD4set))
rc[grep("M1_",rownames(CD4set))]="black"
rc[grep("M2_",rownames(CD4set))]="grey"
rc[grep("M3_",rownames(CD4set))]="yellow"
rc[grep("M4_",rownames(CD4set))]="red"
rc[grep("M5_",rownames(CD4set))]="blue"
rc[grep("M6_",rownames(CD4set))]="cyan"
rc[grep("M7_",rownames(CD4set))]="magenta"
rc[grep("M8_",rownames(CD4set))]="orange"
rc[grep("M9_",rownames(CD4set))]="purple"
rc[grep("M10_",rownames(CD4set))]="brown"
rc[grep("M11_",rownames(CD4set))]="green"
rc[grep("M12_",rownames(CD4set))]="pink"


heatmap.2(CD4set,trace="none",density="none",scale="row",
          cexRow=0.5,cexCol=1.5,  RowSideColors = rc,
          col="redgreen", hclustfun=function(x) hclust(x,method="complete"))


data=CD4set #reg#early.activated#nuocyte#liora.data.all
gene="T helper cells"#"UBC"#"RORA"#"UBC"

color = rep("white", 420)
color[which(CD4set[,"EIF2B1"] > -9 )] = "black"   


d = diffMapFlor(data,ndims=3)
plot3d(x=d$vectors[,1], y=d$vectors[,2], z=d$vectors[,3], xlab="PC 1", ylab="PC 2",zlab="PC 3",col=color, main=gene, size=0.5, type="s")



