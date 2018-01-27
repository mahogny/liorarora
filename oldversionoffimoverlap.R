
######################################################################
### Compare with RORA motifs, any particular direct targets? ---- new version
######################################################################


dat <- read.csv("meme/fimo_nds_1000/fimo.txt",sep="\t",stringsAsFactors = FALSE)
dat <- rbind(dat,read.csv("meme/fimo_nus_1000/fimo.txt",sep="\t",stringsAsFactors = FALSE))

dat <- read.csv("meme/fimo_nds_2000/fimo.txt",sep="\t",stringsAsFactors = FALSE)
dat <- rbind(dat,read.csv("meme/fimo_nus_2000/fimo.txt",sep="\t",stringsAsFactors = FALSE))


genesWithMotif <- dat$sequence_name[dat$X..motif_id=="MA0072.1" & dat$p.value<1e-4]
genesWithMotifSym <- ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% genesWithMotif]

genesWithMotif <- dat$sequence_name[dat$X..motif_id=="MA0071.1" & dat$p.value<1e-4]
genesWithMotifSym <- ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% genesWithMotif]

#intersect(rownames(de_rora)[de_rora$p_val<0.1], genesWithMotifSym)
intersect(rownames(de_rora)[de_rora$p_val<0.01], genesWithMotifSym)






######################################################################
### Compare with RORA motifs, any particular direct targets? ---- older version
######################################################################


#genesWithMotif <- na.omit(rownames(out)[out$m1>=1])
genesWithMotif <- na.omit(rownames(out)[out$m2>=1])  #stricter
genesWithMotifSym <- ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% genesWithMotif]

intersect(rownames(de_rora)[de_rora$p_val<0.1], genesWithMotifSym)

##########################
####### Annotation of the rora DE genes with motifs, using motif m2:

#### DE p-val 0.01
#Slc25a19  higher in Th2 with interesting time dynamics. high early and late
#the human mitochondrial thiamin pyrophosphate transporter (hMTPPT; SLC25A19)
#many papers, metabolic obviously

#S100a4  high in naive, disappears mature. lower Th2. S100 calcium binding protein A4 .
#mouse has immune phenotypes http://www.informatics.jax.org/marker/MGI:1330282

#### DE p-val < 0.1

#Mcm2 peaks at 24h in Th2.   https://en.wikipedia.org/wiki/MCM2   for cell division
#tons of papers. 

#msh2 peaks at 24h in th2/0. https://en.wikipedia.org/wiki/MSH2   dna mismatch repair
#tons of papers

#Snrk lower in Th2. mainly in naive, then going down  https://www.ncbi.nlm.nih.gov/pubmed/?term=snrk 
#https://www.ncbi.nlm.nih.gov/pubmed/23520131 Snrk regulated by lipids. fits role of rora
#Snrk mouse has several metabolic phenotypes http://www.informatics.jax.org/diseasePortal/popup?isPhenotype=true&markerID=MGI:108104&header=homeostasis/metabolism
#http://www.uniprot.org/uniprot/Q8VDU5  snrk claimed to be involved in hematopoiesis and neurons. these are also tissues rora is involved in

##########################
####### Annotation of the rora DE genes with motifs, using motif m1:
### Generally different genes but this motif is less precise. S100a4 agrees. more genes come out

