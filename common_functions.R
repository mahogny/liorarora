library(stringr)
library(topGO)
library(GO.db)
library(org.Mm.eg.db)
library(biomaRt)
library(sqldf)


normalizesym <- function(s) paste(str_sub(s,1,1),str_to_lower(str_sub(s,2)),sep="")


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





#######################################################################
################## Common GO functions ################################
#######################################################################


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
#colnames(human_ensconvert) <- c("ensembl")




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
ortho_mouse_human <- read.csv("mouse_human_ortholog.csv",stringsAsFactors = FALSE)
colnames(ortho_mouse_human) <- c("ens_mouse","ens_human")
ortho_mouse_human <- ortho_mouse_human[ortho_mouse_human$ens_human!="" & ortho_mouse_human$ens_mouse!="",]
## Only 1-1 mappings
ortho_mouse_human_unique <- ortho_mouse_human[
  isUnique(ortho_mouse_human$ens_human) &
    isUnique(ortho_mouse_human$ens_mouse),]
sum(duplicated(ortho_mouse_human_unique$ens_mouse))
sum(duplicated(ortho_mouse_human_unique$ens_human))
ortho_mouse_human




##########################################################################
############### no fan of below #########################################
##########################################################################





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

human_ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart=human_ensembl),stringsAsFactors=FALSE)
colnames(human_ensconvert)[2] <- "mgi_symbol"  #should use a different name. genesymbol
v<-sqldf("select *,count(mgi_symbol) as c from human_ensconvert group by ensembl_gene_id")
human_ensconvert <- human_ensconvert[human_ensconvert$ensembl_gene_id %in% v$ensembl_gene_id[v$c==1],] #there are some bastards including CCL3L3. this one better inserted b$
human_ensconvert$mgi_symbol <- normalizesym(human_ensconvert$mgi_symbol)

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

nametogenesym<-function(x){
  names(x)<-togenesym(names(x))
  x
}

mtogenesym<-function(x){
  colnames(x)<-togenesym(colnames(x))
  rownames(x)<-togenesym(rownames(x))
  x
}



######################################################################
### Read ensembl - transcripts #####################################################
######################################################################
ens_transc_gene <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ens_transc_gene) <- ens_transc_gene$ensembl_transcript_id



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


########################################
## Orthology table human<->mouse
ortho_mouse_human <- read.csv("mouse_human_ortholog.csv",stringsAsFactors = FALSE)
colnames(ortho_mouse_human) <- c("ens_mouse","ens_human")
ortho_mouse_human <- ortho_mouse_human[ortho_mouse_human$ens_human!="" & ortho_mouse_human$ens_mouse!="",]
## Only 1-1 mappings  #good enough for us!
ortho_mouse_human_unique <- ortho_mouse_human[
  isUnique(ortho_mouse_human$ens_human) &
    isUnique(ortho_mouse_human$ens_mouse),]
sum(duplicated(ortho_mouse_human_unique$ens_mouse))
sum(duplicated(ortho_mouse_human_unique$ens_human))


###############################
## Normalization by deseq size factors
normalizeDeseqByGenes <- function(cnt) t(t(cnt)/estimateSizeFactorsForMatrix(cnt[-grep("ERC",rownames(cnt)),]))
normalizeDeseqByERC   <- function(cnt) t(t(cnt)/estimateSizeFactorsForMatrix(cnt[grep("ERC",rownames(cnt)),]))
normalizeDeseqByAll   <- function(cnt) t(t(cnt)/estimateSizeFactorsForMatrix(cnt))
