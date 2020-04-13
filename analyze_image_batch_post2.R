library(stringr) 
library(ggplot2)
library(sqldf)

#####################################
## Make a ggplot multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#######################################################################
################## Read mouse annotations #############################
#######################################################################

#experiment 1 is 30 points
#experiment 2 is 29 points
#wtf?

mousecond <- read.csv("imageanalysis/allslidemeta.csv", stringsAsFactors = FALSE)
mousecond$name <- sprintf("%s-%s", mousecond$expnumber,mousecond$slidename)
mousecond$name <- str_to_lower(mousecond$name)

#mousecond$genotype[mousecond$genotype=="+/+" & mousecond$iscre=="YES"] <- "cre +/+"
#mousecond$genotype[mousecond$genotype=="+/+" & mousecond$iscre=="NO" ] <- "WT"

#######################################################################
################## Read histograms ####################################
#######################################################################

interph <- NULL
allhist <- list()
#nholes <- list()
for(fname in list.files("imageanalysis/images.out",pattern = "*emph*")){
  fs <- str_split_fixed(fname,"_",2)
  fs2 <- str_split_fixed(fs[2],pattern = "\\.",2)[1]
  fs1 <- fs[1]
  #fs[1]
  fs1
  fs2
  
  s <- read.table(sprintf("imageanalysis/images.out/%s",fname))[,1]
  nholes <- length(s)
  s <- sort(s)
  # s <- s[-c(1:10000)]
  # s <- rev(s)[-c(1:50)]
#  s <- s[s>300]
  n <- sprintf("%s-%s",fs1,fs2)
  allhist[[n]] <- s
  
  #Filter anything?
  allperc <- (0:100)/100
  
  onei <- data.frame(
    name=fs1, 
    num=fs2, 
    t(quantile(s,allperc)), 
    nholes=nholes,
    stringsAsFactors = FALSE)
  interph <- rbind(interph,onei)
}


interph

#######################################################################
################## Read overall meta ##################################
#######################################################################


imagemeta <- read.csv("imageanalysis/meta.out.csv", stringsAsFactors = FALSE, header = FALSE, sep=" ")
colnames(imagemeta) <- c("name","dir","sum_nuc","sum_cyto")

#str_split_fixed(imagemeta)

  

#######################################################################
################## Stats - holesize variance, kind of #################
#######################################################################
############ Plot hole size for certain interval vs another interval. deal with compression artifacts?



dat <- merge(interph, mousecond[,c("name","genotype","group","iscre","infected")])
dat$genotype <- factor(dat$genotype, levels = c("+/+","cre +/+","+/flox","flox/flox"))
dat$repr_hole <- dat$X99. / dat$X80.   

#dat <- sqldf("select name, genotype, `group`, avg(repr_hole) as repr_hole, iscre, infected from dat group by name")


make_rel_holesize_plot <- function(group, title=group, iscre="YES", infected="yes"){
  subdat <- dat[dat$group %in% group & dat$iscre==iscre & dat$infected==infected,]
  ggp <- ggplot(subdat, aes(genotype, repr_hole)) +
    ylim(10,90) + 
    geom_point(size=3) +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle=90, hjust=1)) +
    labs(x = title, y = "Emphysemia")
  ggp  
}

panelw=2

pdf("imageanalysis/out_relhole_infected_cre.pdf", w=4*panelw,h=3)
multiplot(cols=4,
  make_rel_holesize_plot("a","Replicate #1"),
  make_rel_holesize_plot("b","Replicate #2"),
  make_rel_holesize_plot("c","Replicate #3"),
  make_rel_holesize_plot(c("a","b","c"),"#1 + #2 + #3"))
dev.off()

##### What if not infected?
pdf("imageanalysis/out_relhole_notinfected_cre.pdf", w=4*panelw,h=3)
multiplot(cols=4,
          make_rel_holesize_plot("a","Replicate #1", infected = "no"),
          make_rel_holesize_plot("c","Replicate #3", infected = "no"))
dev.off()

##### What if not cre?
pdf("imageanalysis/out_relhole_infected_notcre.pdf", w=4*panelw,h=3)
multiplot(cols=4,
          make_rel_holesize_plot("a","Replicate #1", iscre = "NO"),
          make_rel_holesize_plot(c("a","b","c"),"#1 + #2 + #3", iscre = "NO"))
dev.off()



####### Simple test, all together, +/+ vs flox/flox
t.test(
  dat$repr_hole[dat$genotype=="+/+" & dat$iscre=="YES" & dat$infected=="yes"],
  dat$repr_hole[dat$genotype=="flox/flox" & dat$iscre=="YES"  & dat$infected=="yes"])$p.value
# p=0.16

###### Linear model, the one used
dat$expecteffect <- 0
dat$expecteffect[which(dat$genotype=="+/flox" & dat$iscre=="YES"  & dat$infected=="yes")] <- 1
dat$expecteffect[which(dat$genotype=="flox/flox" & dat$iscre=="YES"  & dat$infected=="yes")] <- 2
themod <- lm(repr_hole~expecteffect+group, dat[dat$infected=="yes",])
anova(themod)
# p=0.018188   if only comparing infected ones - I like this test the most


###### Linear model, not making use of increased # alleles giving stronger effect
dat$expecteffect <- 0
dat$expecteffect[which(dat$genotype=="+/flox" & dat$iscre=="YES"  & dat$infected=="yes")] <- 1
dat$expecteffect[which(dat$genotype=="flox/flox" & dat$iscre=="YES"  & dat$infected=="yes")] <- 1
themod <- lm(repr_hole~expecteffect+group, dat[dat$infected=="yes",])
anova(themod)
# p=0.17  but this does not make use of the trend we see

##Pull out example images of high and low
dat[which.min(dat$repr_hole),] #+/+
dat[which.max(dat$repr_hole),] #flox/flox




 
###################
###################
###################
###################
###################



# 
# ############ Plot #holes/total area. Should be lower for flox/flox --------- meah.
# dat <- merge(imagemeta, interph[,c("name","condition","nholes","include")])
# dat <- dat[dat$include=="YES",]
# dat$rel_holes <- dat$nholes/(dat$sum_nuc+dat$sum_cyto)
# ggp <- ggplot(dat, aes(condition, rel_holes)) +
#   geom_point(size=3) +
#   theme(text = element_text(size=20),
#         axis.text.x = element_text(angle=90, hjust=1)) +
#   labs(x = "Condition", y = "pheno")
# print(ggp)
# 
# 
# ############ Plot relative nuclei/cytoplasm. Might be higher in flox/flox
# dat <- merge(imagemeta, interph[,c("name","condition","nholes","include")])
# dat <- dat[dat$include=="YES",]
# dat$rel_nuclei <- dat$sum_nuc/dat$sum_cyto
# ggp <- ggplot(dat, aes(condition, rel_nuclei)) +
#   geom_point(size=3) +
#   theme(text = element_text(size=20),
#         axis.text.x = element_text(angle=90, hjust=1)) +
#   labs(x = "Condition", y = "pheno")
# print(ggp)
# 
# 
# 
# 
# 
# 
# t.test(
#   interph$X80.[interph$include=="YES" & interph$color=="red"],
#   interph$X80.[interph$include=="YES" & interph$color=="blue"])
# 
# sd(interph$X80.[interph$include=="YES" & interph$color=="red"])
# sd(interph$X80.[interph$include=="YES" & interph$color=="blue"])
# 
# 
# #should really average together replicates here!!!
# 
# wstat <- data.frame(
#   name=interph$name,
#   num=interph$num,
#   stat=interph$X80)
# wstat <- merge(wstat, mousecond)
# 
# 
# wstat$pheno<-"ANY"
# wstat$pheno[wstat$genotype=="flox/flox"] <- "KO"
# wstat$pheno[wstat$genotype=="+/flox" & wstat$group=="a"] <- "ctrl"
# wstat$pheno[wstat$genotype=="+/+"    & wstat$group=="b"] <- "ctrl"
# 
# wstat$group <- factor(wstat$group)
# wstat$pheno <- factor(wstat$pheno)
# 
# themod <- lm(stat~pheno+group, wstat)
# anova(themod)
# 
# 
# 

#######################################################################
################## Stats - holesize ###################################
#######################################################################
############ Plot hole size for certain interval

dat <- merge(interph, mousecond[,c("name","genotype","group","iscre")])
dat$genotype <- factor(dat$genotype, levels = c("+/+","cre +/+","+/flox","flox/flox"))
dat$abs_hole <- dat$X98.

dat <- sqldf("select name, genotype, `group`, avg(abs_hole) as abs_hole, iscre from dat group by name")


make_holesize_plot <- function(group, title=group){
  subdat <- dat[dat$group %in% group & dat$iscre=="YES",]
  ggp <- ggplot(subdat, aes(genotype, abs_hole)) +
    geom_point(size=3) +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle=90, hjust=1)) +
    labs(x = title, y = "Emphysemia")
  ggp  
}

multiplot(cols=2,
          make_holesize_plot("a"),
          make_holesize_plot("b"),
          make_holesize_plot("c"),
          make_holesize_plot(c("a","b","c"),"Combined"))



thegroup <- "a"
t.test(
  dat$abs_hole[dat$group==thegroup & dat$genotype=="WT"],
  dat$abs_hole[dat$group==thegroup & dat$genotype=="flox/flox"])$p.value

thegroup <- "b"
t.test(
  dat$abs_hole[dat$group==thegroup & dat$genotype=="WT"],
  dat$abs_hole[dat$group==thegroup & dat$genotype=="flox/flox"])$p.value

thegroup <- "c"
t.test(
  dat$abs_hole[dat$group==thegroup & dat$genotype=="+/flox"],
  dat$abs_hole[dat$group==thegroup & dat$genotype=="flox/flox"])$p.value

t.test(
  dat$abs_hole[dat$genotype=="+/flox"],
  dat$abs_hole[dat$genotype=="flox/flox"])$p.value
