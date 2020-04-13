library(stringr)
library(sqldf)

#todo. actually correct the layout!
cytolay_f <- read.table("/home/mahogny/Dropbox/applyPI/liora/upstream/layout_cytoscreen_123_corrected.csv", 
                        stringsAsFactors = FALSE, sep = ",")
cytolay_f <- cytolay_f[,c(12,11,10,9,8,7,6,5,4,3,2,1)] #correct for pipeting error...

cytolay <- NULL
for(i in 1:12){
  for(j in 1:8){
    cytolay <- rbind(cytolay,
                     data.frame(well=sprintf("%s%s",c("A","B","C","D","E","F","G","H")[j],i), 
                                treatment=cytolay_f[j,i], stringsAsFactors = FALSE))
  }
}

#TODO also include plate 4


cytolay

dat <- read.csv("/home/mahogny/Dropbox/applyPI/liora/upstream/facsdata.csv", stringsAsFactors = FALSE)
colnames(dat) <- c("Dataset","activated","sc","lowdiv","mouse")
dat$well <- str_split_fixed(dat$Dataset, "_", 5)[,3]
dat$plate <- as.double(str_split_fixed(dat$Dataset, "_", 5)[,2])
dat$mouse[is.na(dat$mouse)] <- dat$plate[is.na(dat$mouse)]

alldat <- merge(dat, cytolay)


comp <- NULL
for(g in unique(alldat$treatment)){
  
  p.lowdiv <- t.test(
    alldat$lowdiv[alldat$treatment=="wt"],
    alldat$lowdiv[alldat$treatment==g]
  )$p.val
  
  fc.lowdiv <- mean(alldat$lowdiv[alldat$treatment==g]) / mean(alldat$lowdiv[alldat$treatment=="wt"])
  
  comp <- rbind(comp,data.frame(
    gene=g,
    p.lowdiv=p.lowdiv, 
    fc.lowdiv=fc.lowdiv))
}

################################################
################################################
################################################  
comp <- comp[order(comp$p.lowdiv),]
comp
#plot(comp$fc.lowdiv)

#but this means il4, il13 and lowdiv all decrease lowdiv!

plot(comp$fc.lowdiv, -comp$p.lowdiv, cex=0)
text(comp$fc.lowdiv, -comp$p.lowdiv, labels = comp$gene)



################################################
################################################  just show it all!
################################################  


pertreat <- sqldf("select well, treatment, avg(lowdiv) as lowdiv from alldat group by well")
#pertreat <- pertreat[order(pertreat$lowdiv),]
pertreat$color <- "black"
pertreat$color[pertreat$treatment=="wt"] <- "blue"

plot(pertreat$lowdiv,cex=0)
text(pertreat$lowdiv, labels = pertreat$treatment, col=pertreat$color)
  
