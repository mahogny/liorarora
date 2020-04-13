#########################################################################
######################## Plot rora over time ############################
#########################################################################

#could supplement with the naive in vitro(?) - what was sorted?

keep <- which(cellcond37$group=="tc")
cellcond37[keep,]

ensid_rora <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Rora"]
enstrans_rora <- ens_transc_gene$ensembl_transcript_id[ens_transc_gene$ensembl_gene_id==ensid_rora]

tc_rora <- data.frame(
  roracount=as.double(cnt37[ensid_rora,keep]),  #sum appropriate, not mean
  day=as.double(cellcond37[keep,]$day))

pdf("out/roratc.pdf",width = 5, height = 3)
plot(tc_rora$day, tc_rora$roracount,
     ylim=c(0,300),pch=19)
tc_rora_mean <- sqldf("select distinct avg(roracount) as c, day from tc_rora group by day")
lines(x=tc_rora_mean$day, y=tc_rora_mean$c,col="red")
dev.off()

#########################################################################
###################### Make a combined table ############################
#########################################################################

################### bulk data
de_merge <- merge(de_overexp, de_tissue, all=TRUE)
de_merge <- merge(de_merge, de_nipo30, all=TRUE)
de_merge <- merge(de_merge, de_nipo0,  all=TRUE)
de_merge <- merge(de_merge, de_treg_spl_0,  all=TRUE)
de_merge <- merge(de_merge, de_nipo_7,  all=TRUE)
de_merge <- merge(de_merge, de_tissue_combined,  all=TRUE)
de_merge <- merge(de_merge, de_th17,  all=TRUE)
de_merge <- merge(de_merge, de_livermacro,  all=TRUE)
de_merge <- merge(de_merge, de_thp1,  all=TRUE)

de_merge <- merge(de_merge, de_skin_treg_ens,  all=TRUE) ##### this adds mgi_symbol


################## single cell data, nipo
allde_sc<-NULL
for(g in c("Treg","Th2","Naive","Th1")){ #"all",
  print(g)
  ### Straight scRNAseq comparison
  onemore <- data.frame(
    fc=de_sc[[g]]$avg_logFC,#log2FoldChange, 
    pvalue=de_sc[[g]]$p_val,#pvalue, 
    mgi_symbol = rownames(de_sc[[g]]),stringsAsFactors = FALSE)
  colnames(onemore) <- c(
    sprintf("fc_sc_%s",g), 
    sprintf("p_sc_%s",g), 
    "mgi_symbol")
  if(is.null(allde_sc)){
    allde_sc <- onemore
  } else {
    allde_sc <- merge(allde_sc, onemore, all=TRUE)
  }
  
  ### in vitro Bulk scRNAseq
  # onemore <- data.frame(
  #   fc=de_bulk[[g]]$log2FoldChange,
  #   pvalue=de_bulk[[g]]$pvalue,
  #   ensembl_gene_id = rownames(de_bulk[[g]]),stringsAsFactors = FALSE)
  # colnames(onemore) <- c(
  #   sprintf("fc_scb_%s",g), 
  #   sprintf("p_scb_%s",g), 
  #   "mgi_symbol")
  # allde_sc <- merge(allde_sc, onemore, all=TRUE)
}
de_merge <- merge(de_merge, allde_sc, all=TRUE)


list_fc_cond <- str_split_fixed(colnames(de_merge)[grep("fc_",colnames(de_merge))],"_",2)[,2]


######## Fill in empty values
for(i in grep("fc_", colnames(de_merge))){
  de_merge[is.na(de_merge[,i]),i] <- 0
}
for(i in grep("p_", colnames(de_merge))){
  de_merge[is.na(de_merge[,i]),i] <- 1
}

######## Fix missing gene symbols. Introduced by previous stupid dataset
for(i in which(is.na(de_merge$mgi_symbol))){
  de_merge$mgi_symbol[i] <- ensconvert$mgi_symbol[ensconvert$ensembl_gene_id==de_merge$ensembl_gene_id[i]][1]
}



#########################################################################
###################### Plot FC comparisons ##############################
#########################################################################


plotFCs <- function(seta,setb, pcut.a=0.1, pcut.b=0.1, cex=0.5, fccut=0.01, outfile="",cexp=0){
  subde <- data.frame(
    mgi_symbol=de_merge$mgi_symbol,
    ensembl_gene_id=de_merge$ensembl_gene_id,
    fc_a=de_merge[,sprintf("fc_%s",seta)],
    fc_b=de_merge[,sprintf("fc_%s",setb)],
    p_a =de_merge[,sprintf("p_%s",seta)],
    p_b =de_merge[,sprintf("p_%s",setb)],
    stringsAsFactors = FALSE)
  
  subde$mgi_symbol[which(is.na(subde$mgi_symbol))] <- subde$ensembl_gene_id[which(is.na(subde$mgi_symbol))]
  
  keep <- abs(subde$fc_a)>fccut & abs(subde$fc_b)>fccut & subde$p_a<pcut.a & subde$p_b<pcut.b
#  print(sum(keep))
 # print(subde[keep,])
  
  
  thec <- cor(
    subde$fc_a[keep],
    subde$fc_b[keep], method = "spearman")
  
#  if(sum(keep)>20 & thec>0.2){
    if(outfile!=""){
      if(endsWith(outfile,"pdf")){
        pdf(outfile)
      }
      if(endsWith(outfile,"png"))
        png(outfile, width=1600, height = 1200)
    }
    plot(subde$fc_a[keep], subde$fc_b[keep], 
         cex=cexp, pch=19, xlab=seta, ylab=setb)
    lines(c(-30,30),c(0,0), col="gray")
    lines(c(0,0),c(-30,30), col="gray")
    text(subde$fc_a[keep], subde$fc_b[keep], 
         label=subde$mgi_symbol[keep], cex=cex)
    if(outfile!=""){
      dev.off()
    }
  #}  
}


#could filter out lowly exp genes here. or color by exp

write.csv(de_merge,"out/merged_all_de.csv")


for(i in list_fc_cond){
  for(j in list_fc_cond){
    if(i!=j){
#      i <- "colon"
 #     j <- "lung"
      
      #not done yet: should only keep unique ones
      
      plotFCs(i,j, outfile=sprintf("out.de_comp/allvsall_%s_%s.pdf",i,j),cex=1, cexp=0.5) 
#      plotFCs(i,j, outfile=sprintf("out.de_comp/allvsall_%s_%s.png",i,j),cex=2) 
    }
  }
}


#plotFCs(i,j, outfile=sprintf("out.de_comp/allvsall_%s_%s.pdf",i,j),cex=1, cexp=0.5) 


#########################################################################
###################### Confirm DE genes on the SC level #################
#########################################################################


de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol=="S100a4",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg"))]
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol=="S100a4",sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg"))]


de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol=="Arntl",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg"))]
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol=="Arntl",sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg"))]


theg <- "Cxcr6"
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg"))]
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg"))]

theg <- "Il10"
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg"))]
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg"))]


theg <- "Nebl"
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg"))]
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg"))]

theg <- "Tox2"
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg"))]
de_merge[!is.na(de_merge$mgi_symbol) & de_merge$mgi_symbol==theg,sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg"))]


##naive
desummary_p <- de_merge[
  de_merge$mgi_symbol %in% keepgenes & de_merge$fc_overexp*de_merge$fc_sc_Naive>0 & de_merge$p_sc_Naive<0.01,
  c("mgi_symbol","p_overexp",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p
###=>   Chd3   Gimap7   Jak3   Ssb?

##Th1
desummary_p <- de_merge[
  de_merge$mgi_symbol %in% keepgenes & de_merge$fc_overexp*de_merge$fc_sc_Th1>0 & de_merge$p_sc_Th1<0.01,
  c("mgi_symbol","p_overexp",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p
###=> Mcm6?

##Th2
desummary_p <- de_merge[
  de_merge$mgi_symbol %in% keepgenes & de_merge$fc_overexp*de_merge$fc_sc_Th2>0 & de_merge$p_sc_Th2<0.01,
  c("mgi_symbol","p_overexp",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p
###=> Npm1  Ppm1d
# mgi_symbol p_overexp p_sc_Naive  p_sc_Th1     p_sc_Th2 p_sc_Treg
# 2333       Anxa6 0.3872788 0.33368119 0.3511690 2.849297e-03 1.0000000
# 3465       Bscl2 0.4376718 1.00000000 0.6911929 3.198641e-03 1.0000000
# 5425      Ctdsp1 0.6004154 1.00000000 0.7445856 3.238539e-03 1.0000000
# 15241     Mrpl43 0.7759774 1.00000000 1.0000000 3.704113e-03 1.0000000
# 15630       Nab1 0.4338973 0.04740998 1.0000000 4.998633e-03 0.2662996
# 15674       Nars 0.8712107 1.00000000 0.5491130 3.978811e-03 1.0000000
# 16144       Npm1 0.4458613 1.00000000 0.5818750 9.119007e-05 0.2140565 ***
# 19094      Ppm1d 0.1667709 1.00000000 1.0000000 1.003628e-03 1.0000000 ***
# 20803      Rplp2 0.5680439 1.00000000 0.7808477 3.638145e-03 0.3080076
# 23697     Tgoln1 0.9948183 0.74423046 1.0000000 2.711756e-03 1.0000000
# 26636      Zfp53 0.5900930 1.00000000 0.3860046 1.471582e-03 1.0000000

#de_merge$p_overexp
##Treg
desummary_p <- de_merge[
  de_merge$mgi_symbol %in% keepgenes & de_merge$fc_overexp*de_merge$fc_sc_Treg>0 & de_merge$p_sc_Treg<0.01,
  c("mgi_symbol","p_overexp","p_skin_treg",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p

# mgi_symbol  p_overexp p_sc_Naive  p_sc_Th1   p_sc_Th2    p_sc_Treg
# 2611     Arl6ip5 0.13452734  0.3771406 1.0000000 1.00000000 0.0032033572
# 3784      Camk2b 0.90610402  1.0000000 1.0000000 1.00000000 0.0019535985
# 4294       Cdc37 0.83224586  0.4982072 1.0000000 1.00000000 0.0022131172
# 4664        Chd6 0.92951139  1.0000000 0.5461574 1.00000000 0.0044496600
# 7169      Exosc1 0.83770137  1.0000000 1.0000000 1.00000000 0.0002735498
# 11582      Ift27 0.69724480  1.0000000 1.0000000 1.00000000 0.0015445672
# 12192      Kdm4c 0.05932195  0.7245137 0.9360020 1.00000000 0.0006074434 ******
# 15632      Nabp1 0.07274987  0.2276460 0.9290726 0.07438566 0.0037227776 *****
# 19725      Psma2 0.43425834  0.5418203 1.0000000 1.00000000 0.0012251279
# 21659      Skap1 0.35960847  1.0000000 1.0000000 1.00000000 0.0008707014
# 22509      Sorl1 0.07009682  1.0000000 1.0000000 1.00000000 0.0004667415 **** 
# 22811      Srsf2 0.78192916  0.3062243 1.0000000 1.00000000 0.0029760850
# 24984       Ubn2 0.70832395  0.9036861 0.5002101 1.00000000 0.0018821680
# 25380      Usp40 0.95296682  1.0000000 1.0000000 1.00000000 0.0012898600


##Treg
desummary_p <- de_merge[
  de_merge$mgi_symbol %in% keepgenes & 
    de_merge$fc_overexp*de_merge$fc_sc_Treg>0 & 
    de_merge$fc_skin_treg*de_merge$fc_sc_Treg>0 & de_merge$p_sc_Treg<0.01,
  c("mgi_symbol","p_overexp","p_skin_treg",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p



#########################################################################
###################### Heatmap of FC ####################################
#########################################################################

# 
# sum(keepgene)
# 
# 
# keepgenes <- c(
#   de_merge$mgi_symbol[de_merge$p_sc_Naive<0.05 & abs(de_merge$fc_sc_Naive)>=1],
#   de_merge$mgi_symbol[de_merge$p_sc_Th1<0.05   & abs(de_merge$fc_sc_Th1)>=1],
#   de_merge$mgi_symbol[de_merge$p_sc_Th2<0.05   & abs(de_merge$fc_sc_Th2)>=1],
#   de_merge$mgi_symbol[de_merge$p_sc_Treg<0.05  & abs(de_merge$fc_sc_Treg)>=1])
# keepgenes
# #No overlap!
# 
# 
# fccut=0.3
# pcut=0.002
# keepgenes <- c(
#   de_merge$mgi_symbol[de_merge$p_sc_Naive<pcut & abs(de_merge$fc_sc_Naive)>=fccut],
#   de_merge$mgi_symbol[de_merge$p_sc_Th1<pcut   & abs(de_merge$fc_sc_Th1)>=fccut],
#   de_merge$mgi_symbol[de_merge$p_sc_Th2<pcut   & abs(de_merge$fc_sc_Th2)>=fccut],
#   de_merge$mgi_symbol[de_merge$p_sc_Treg<0.1   & abs(de_merge$fc_sc_Treg)>=0.1])
# keepgenes


fccut=0.5
pcut=0.0001
keepgenes <- sort(unique(c(
  de_merge$mgi_symbol[de_merge$p_sc_Naive<0.0001   & abs(de_merge$fc_sc_Naive)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Th1  <0.00001  & abs(de_merge$fc_sc_Th1)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Th2  <0.001    & abs(de_merge$fc_sc_Th2)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Treg <0.0001   & abs(de_merge$fc_sc_Treg)>=fccut])))
sort(keepgenes)

fccut=0
pcut=0.0001
keepgenes <- sort(unique(c(
  de_merge$mgi_symbol[de_merge$p_sc_Naive<0.001   & abs(de_merge$fc_sc_Naive)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Th1  <0.001  & abs(de_merge$fc_sc_Th1)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Th2  <0.001    & abs(de_merge$fc_sc_Th2)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Treg <0.001   & abs(de_merge$fc_sc_Treg)>=fccut])))
sort(keepgenes)


fccut=0
pcut=0.0001
keepgenes <- sort(unique(c(
  de_merge$mgi_symbol[de_merge$p_sc_Naive<0.005   & abs(de_merge$fc_sc_Naive)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Th1  <0.005  & abs(de_merge$fc_sc_Th1)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Th2  <0.005    & abs(de_merge$fc_sc_Th2)>=fccut],
  de_merge$mgi_symbol[de_merge$p_sc_Treg <0.005   & abs(de_merge$fc_sc_Treg)>=fccut])))
sort(keepgenes)



#de_merge[]
#de_merge$fc_overexp*de_merge$fc_sc_Naive>0 & de_merge$p_sc_Naive<0.01

desummary_p <- de_merge[
  de_merge$mgi_symbol %in% keepgenes & de_merge$fc_overexp*de_merge$fc_sc_Naive>0 & de_merge$p_sc_Naive<0.01,
  c("mgi_symbol",sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p

#de_merge$mgi_symbol


###can use side-color to denote which cell type they are DE in

desummary_p <- de_merge[de_merge$mgi_symbol %in% keepgenes,c(sprintf("p_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_p<0.01

desummary_fc <- de_merge[de_merge$mgi_symbol %in% keepgenes,c(sprintf("fc_sc_%s",c("Naive","Th1","Th2","Treg")))]
desummary_fc$mgi_symbol <- de_merge$mgi_symbol[de_merge$mgi_symbol %in% keepgenes]

desummary_fc <- unique.array(desummary_fc)
rownames(desummary_fc) <- desummary_fc$mgi_symbol
desummary_fc <- desummary_fc[,-5]
colnames(desummary_fc) <- c("Naive","Th1","Th2","Treg")

pdf("out.niposs2/new_topDE_heatmapFC_2.pdf", width = 3, height = 6)
heatmap.2(
  #symm = TRUE,
  col=colorRampPalette(c("red", "white", "green"))(n = 299),
  cexRow = 1.5,
  as.matrix(desummary_fc),
  density.info="none",
  trace="none")
dev.off()




de_merge$mgi_symbol[de_merge$p_sc_Naive<pcut & abs(de_merge$fc_sc_Naive)>=fccut]


#########################################################################
###################### Compare with FIMO ################################
#########################################################################


#########################################################################
###################### GO analysis ######################################
#########################################################################

stopgosym(
  unique(togenesym(de_merge$ensembl_gene_id[rank(de_merge$p_sc_Treg)<200])),
  unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_merge$ensembl_gene_id]))
# 
# GO.ID                                   term.name Annotated Significant Expected p.value
# 1  GO:0007600                          sensory perception      1939        1933  1924.50   0.008
# 2  GO:0006082              organic acid metabolic process       917         915   910.14   0.030
# 3  GO:0043436                   oxoacid metabolic process       904         902   897.24   0.032
# 4  GO:0007606     sensory perception of chemical stimulus      1423        1418  1412.36   0.041
# 5  GO:0019752           carboxylic acid metabolic process       844         842   837.69   0.046
# 6  GO:0050865               regulation of cell activation       638         637   633.23   0.046
# 7  GO:0022407            regulation of cell-cell adhesion       392         392   389.07   0.051
# 8  GO:0050673               epithelial cell proliferation       385         385   382.12   0.054
# 9  GO:0002694          regulation of leukocyte activation       598         597   593.53   0.059
# 10 GO:0050678 regulation of epithelial cell proliferat...       323         323   320.58   0.087
# 11 GO:0051249         regulation of lymphocyte activation       522         521   518.10   0.095
# 12 GO:0007186 G-protein coupled receptor signaling pat...      1870        1861  1856.01   0.099
# 13 GO:1903037 regulation of leukocyte cell-cell adhesi...       301         301   298.75   0.103
# 14 GO:0006909                                phagocytosis       293         293   290.81   0.109
# 15 GO:0070371                       ERK1 and ERK2 cascade       293         293   290.81   0.109
# 16 GO:0050863             regulation of T cell activation       286         286   283.86   0.115
# 


stopgosym(
  unique(togenesym(de_merge$ensembl_gene_id[rank(de_merge$p_sc_Naive)<200])),
  unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_merge$ensembl_gene_id]))


# GO.ID                                   term.name Annotated Significant Expected p.value
# 1  GO:0061024                       membrane organization       735         734   729.70   0.029
# 2  GO:0002253               activation of immune response       443         443   439.81   0.039
# 3  GO:0001525                                angiogenesis       443         443   439.81   0.039
# 4  GO:0007610                                    behavior       675         674   670.13   0.042
# 5  GO:0001944                     vasculature development       673         672   668.15   0.043
# 6  GO:1901615 organic hydroxy compound metabolic proce...       424         424   420.94   0.045
# 7  GO:0001568                    blood vessel development       643         642   638.36   0.052
# 8  GO:0002764 immune response-regulating signaling pat...       396         396   393.14   0.056
# 9  GO:0002757 immune response-activating signal transd...       383         383   380.24   0.061
# 10 GO:0051090 regulation of DNA binding transcription ...       358         358   355.42   0.073
# 11 GO:0006631                fatty acid metabolic process       354         354   351.45   0.076
# 12 GO:0002237    response to molecule of bacterial origin       329         329   326.63   0.091
# 13 GO:0048514                  blood vessel morphogenesis       546         545   542.06   0.093
# 14 GO:0032787       monocarboxylic acid metabolic process       526         525   522.21   0.104
# 15 GO:0032496              response to lipopolysaccharide       310         310   307.76   0.104
# 16 GO:0050890                                   cognition       292         292   289.89   0.119



stopgosym(
  unique(togenesym(de_merge$ensembl_gene_id[rank(de_merge$p_sc_Th2)<200])),
  unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% de_merge$ensembl_gene_id]))

# GO.ID                                   term.name Annotated Significant Expected p.value
# 1  GO:0008285 negative regulation of cell proliferatio...       651         650   646.34    0.05
# 2  GO:0006959                     humoral immune response       292         292   289.91    0.12
# 3  GO:0051170                              nuclear import       282         282   279.98    0.13
# 4  GO:0006606                 protein import into nucleus       278         278   276.01    0.13
# 5  GO:0005996            monosaccharide metabolic process       244         244   242.25    0.17
# 6  GO:0050851 antigen receptor-mediated signaling path...       243         243   241.26    0.17
# 7  GO:0070663       regulation of leukocyte proliferation       227         227   225.37    0.19
# 8  GO:1904950 negative regulation of establishment of ...       225         225   223.39    0.20
# 9  GO:0010506                     regulation of autophagy       223         223   221.40    0.20
# 10 GO:0051224    negative regulation of protein transport       221         221   219.42    0.20
# 11 GO:0051169                           nuclear transport       411         410   408.06    0.20
# 12 GO:0002696 positive regulation of leukocyte activat...       410         409   407.06    0.20
# 13 GO:0060291             long-term synaptic potentiation       219         219   217.43    0.21
# 14 GO:0032944 regulation of mononuclear cell prolifera...       219         219   217.43    0.21
# 15 GO:0006913                 nucleocytoplasmic transport       408         407   405.08    0.21
# 16 GO:0050670      regulation of lymphocyte proliferation       217         217   215.45    0.21
# 17 GO:0046822 regulation of nucleocytoplasmic transpor...       216         216   214.45    0.21
# 18 GO:0019318                    hexose metabolic process       214         214   212.47    0.21
# 19 GO:0070302 regulation of stress-activated protein k...       208         208   206.51    0.22
# 20 GO:0032872 regulation of stress-activated MAPK casc...       207         207   205.52    0.22
# 21 GO:0050663                          cytokine secretion       202         202   200.55    0.23
# 22 GO:0031668 cellular response to extracellular stimu...       200         200   198.57    0.24
# 23 GO:0034504             protein localization to nucleus       382         381   379.26    0.24
# 24 GO:0042098                        T cell proliferation       197         197   195.59    0.24
# 25 GO:0042742               defense response to bacterium       377         376   374.30    0.24
# 26 GO:0032984          macromolecular complex disassembly       191         191   189.63    0.25
# 27 GO:0030258                          lipid modification       190         190   188.64    0.25
# 28 GO:0051251 positive regulation of lymphocyte activa...       370         369   367.35    0.25
# 29 GO:0097696                                STAT cascade       189         189   187.65    0.26
# 30 GO:0007259                            JAK-STAT cascade       189         189   187.65    0.26
# 31 GO:0072594 establishment of protein localization to...       536         534   532.16    0.26
# 32 GO:0032386       regulation of intracellular transport       534         532   530.17    0.26
# 33 GO:0072376                  protein activation cascade       186         186   184.67    0.26
# 34 GO:0007254                                 JNK cascade       185         185   183.67    0.26
# 35 GO:0006006                   glucose metabolic process       185         185   183.67    0.26
# 36 GO:0006302                  double-strand break repair       184         184   182.68    0.26
# 37 GO:0006974    cellular response to DNA damage stimulus       687         684   682.08    0.27
# 38 GO:0007166     cell surface receptor signaling pathway      2364        2350  2347.06    0.27
# 39 GO:0030307          positive regulation of cell growth       179         179   177.72    0.27
# 40 GO:1904589                regulation of protein import       178         178   176.72    0.28
# 41 GO:0050707            regulation of cytokine secretion       178         178   176.72    0.28
# 42 GO:0043200                      response to amino acid       177         177   175.73    0.28
# 43 GO:0006914                                   autophagy       352         351   349.48    0.28
# 44 GO:0061919      process utilizing autophagic mechanism       352         351   349.48    0.28
# 45 GO:0006836                  neurotransmitter transport       175         175   173.75    0.28
# 46 GO:0002274                myeloid leukocyte activation       175         175   173.75    0.28
# 47 GO:0016236                              macroautophagy       175         175   173.75    0.28
# 48 GO:0042306 regulation of protein import into nucleu...       175         175   173.75    0.28
# 49 GO:0031669        cellular response to nutrient levels       174         174   172.75    0.28
# 50 GO:0006956                       complement activation       173         173   171.76    0.29
# 51 GO:0000075                       cell cycle checkpoint       172         172   170.77    0.29
# 52 GO:1902532 negative regulation of intracellular sig...       507         505   503.37    0.29
# 53 GO:0002455 humoral immune response mediated by circ...       169         169   167.79    0.30
# 54 GO:0046328                   regulation of JNK cascade       169         169   167.79    0.30
# 55 GO:0046425              regulation of JAK-STAT cascade       168         168   166.80    0.30
# 56 GO:1904892                  regulation of STAT cascade       168         168   166.80    0.30
# 57 GO:0022411              cellular component disassembly       338         337   335.58    0.30
# 58 GO:0050853           B cell receptor signaling pathway       166         166   164.81    0.30
# 59 GO:0071230    cellular response to amino acid stimulus       165         165   163.82    0.30
# 60 GO:0007219                     Notch signaling pathway       165         165   163.82    0.30
# > 




#########################################################################
###################### Compare with chip ################################
#########################################################################

dim(anno_rora)  #only 596 genes. be careful with overlap

de_merge_chip <- de_merge
de_merge_chip <- de_merge_chip[de_merge_chip$ensembl_gene_id %in% anno_rora$Nearest.Ensembl,]

de_merge_chip <- de_merge_chip[order(de_merge_chip$p_th17),]
de_merge_chip$mgi_symbol[1:50]


de_merge_chip <- de_merge_chip[order(de_merge_chip$p_spleen),]
de_merge_chip$mgi_symbol[1:50]

de_merge_chip <- de_merge_chip[order(de_merge_chip$p_overexp),]
de_merge_chip$mgi_symbol[1:50]

de_merge_chip <- de_merge_chip[order(de_merge_chip$p_colon),]
de_merge_chip$mgi_symbol[1:50]

de_merge_chip <- de_merge_chip[order(de_merge_chip$p_gut),]
de_merge_chip$mgi_symbol[1:50]

de_merge_chip <- de_merge_chip[order(de_merge_chip$p_macro),]
de_merge_chip$mgi_symbol[1:50]
#Trim24 etc, 

intersect(
  de_merge_chip[order(de_merge_chip$p_spleen),]$mgi_symbol[1:80],
  de_merge_chip[order(de_merge_chip$p_th17),]$mgi_symbol[1:80])
#"Mcph1"  "Trim24"
#"Mcph1"  "Pdcd4"  "Trim24" "Nav2"   "Rap2a" 

#[1] "Rbl2"   "Mcph1"  "Pdcd4"  "Trim24" "Nav2"   "Rap2a"  "Rhbdd1" "Prc1"   "Ccna2"  "Rtca"   "Mthfr"  "Lrrc28" "Rad21" 

## Note: Rbl2 comes up in oe5_rora as well


intersect(
  de_merge_chip[order(de_merge_chip$p_overexp),]$mgi_symbol[1:80],
  de_merge_chip[order(de_merge_chip$p_th17),]$mgi_symbol[1:80])
#[1] "Mcph1"  "Trim24" "Sorl1" 

#"Gcc2"   "Mcph1"  "Trim24" "Sorl1"  "Lnpep"  "Ccna2"

# [1] "Gcc2"    "Mcph1"   "Trim24"  "Rbl2"    "Vcpip1"  "Sorl1"   "Lnpep"   "Rhbdd1"  "Ccna2"   "Incenp"  "Thoc1"   "Rps6ka2" "Txndc16" "Mmd"     "Zfp438"  "Dck"    
# [17] "Pttg1"   "Zfp277" 


intersect(
  de_merge_chip[order(de_merge_chip$p_overexp),]$mgi_symbol[1:80],
  de_merge_chip[order(de_merge_chip$p_macro),]$mgi_symbol[1:80])
#Cmtm6"  "Trim24" "Nampt" 

#[1] "Cmtm6"   "Trim24"  "Nampt"   "Dcun1d4" "Klf11"   "Bcl10"   "Srd5a3"  "Tle4"    "Mis12"   "Mmd"     "Hnrnpk"  "Pttg1"  




# "Mcph1"    slightly DE. expressed and increasing. https://en.wikipedia.org/wiki/Microcephalin    can explain rora and brain
#           http://www.mousephenotype.org/data/genes/MGI:2443308    metabolic and brain phenotype
#           http://www.informatics.jax.org/diseasePortal/genoCluster/view/36376  more auto-immune diabetes. higher sensitivity to radiation.
# "Trim24"  http://www.informatics.jax.org/marker/MGI:109275  can connect to liver cancer (proliferation). 

# Ccna2   Cell cycle control

# Rhbdd1  expressedish. secretion and apoptosis. cleaves proteins in the membrane


write.csv(de_merge_chip,"out/merged_all_de_chip.csv")




