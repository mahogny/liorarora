library(stringr) 
library(ggplot2)

#######################################################################
################## Read mouse annotations #############################
#######################################################################

#experiment 1 is 30 points
#experiment 2 is 29 points
#wtf?

mousecond <- read.csv("imageanalysis/allslidemeta.csv", stringsAsFactors = FALSE)
mousecond$name <- sprintf("%s-%s", mousecond$expnumber,mousecond$slidename)

mousecond$genotype[mousecond$genotype=="+/+" & mousecond$iscre=="YES"] <- "cre +/+"
mousecond$genotype[mousecond$genotype=="+/+" & mousecond$iscre=="NO" ] <- "WT"

#######################################################################
################## Read histograms ####################################
#######################################################################

interph <- NULL
allhist <- list()
#nholes <- list()
for(fname in list.files("images.out",pattern = "*emph*")){
  fs <- str_split_fixed(fname,"_",2)
  fs2 <- str_split_fixed(fs[2],pattern = "\\.",2)[1]
  fs1 <- fs[1]
  #fs[1]
  fs1
  fs2
  
  s <- read.table(sprintf("images.out/%s",fname))[,1]
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


imagemeta <- read.csv("images.meta.csv", stringsAsFactors = FALSE, header = FALSE, sep=" ")
colnames(imagemeta) <- c("name","dir","sum_nuc","sum_cyto")

#str_split_fixed(imagemeta)


#######################################################################
################## Stats ##############################################
############################################################### ########


#color by replicate?

mousecond$name <- str_to_lower(mousecond$name)
#mousecond$name[mousecond$genotype=="flox/flox"]

mousecond$genotype

#+/flox

doA <- TRUE


## Label KO controls
#mousecond$name[mousecond$genotype=="+/flox" & str_sub(interph$name,1,6) =="q31.16"] <- 
interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/flox"]]    <- "blue"

# if(doA){
#   #First group
#   interph$color <- "black"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="flox/flox"]] <- "red"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/flox"]]    <- "blue"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/+"]]       <- "green"
# 
#   interph$condition <- "nocond"
#   interph$condition[interph$color=="red"] <- "flox/flox"
#   interph$condition[interph$color=="blue"] <- "+/flox"
#   interph$condition[interph$color=="green"] <- "+/+"
#   
#   interph$include <- "NO"
#   interph$include[
#     interph$color!="black" & 
#     str_sub(interph$name,1,6) =="q31.16"] <- "YES"
#   
# } else {
#   ########## they call it something else? see plot #2. included this data? or is it part of #1?
#   
#   
#   #First group
#   interph$color <- "black"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="flox/flox"]] <- "red"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/+"]]       <- "blue"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/+"]]       <- "green"
#   
#   # interph$include[interph$color!="black"] <- "YES"
#   
#   interph$include <- "NO"
#   interph$include[
#     interph$color!="black" & 
#       str_sub(interph$name,1,6) =="q31.17"] <- "YES"
# }
# 
# #interph[,c("include","color","name")]
# # 
# # ############ Plot histogram over hole sizes
# # plot(c(1:101),log(1+as.double(interph[1,c(1:101)+2,2])),type="l",ylim=c(0,10),col="white")
# # for(i in 1:nrow(interph)){
# #   if(interph$include[i]=="YES"){
# #     lines(c(1:101),log(1+as.double(interph[i,c(1:101)+2,2])), col=interph$color[i])
# #   }
# # }
# 

#mousecond$

dat <- merge(interph, mousecond[,c("name","genotype","group")])
#interph <- mousecond[,c("name","genotype","group")]

mousecond$name
#mousecond$

# interph$include <- "NO"
# interph$include[interph$name %in% mousecond$name[mousecond$group=="a"]] <- "YES"
# 
# interph$condition <- ""
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="flox/flox"]] <- "red"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/+"]]       <- "blue"
#   interph$color[interph$name %in% mousecond$name[mousecond$genotype=="+/+"]]       <- "green"
#

############ Plot hole size for certain interval
dat <- merge(imagemeta, interph[,c("name","X99.","include")])
dat <- dat[dat$group=="a",]
dat$repr_hole <- dat[,grep("X",colnames(dat))]
ggp <- ggplot(dat, aes(condition, repr_hole)) +
  geom_point(size=3) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "Condition", y = "pheno")
print(ggp)


############ Plot hole size for certain interval vs another interval. deal with compression artifacts?
dat <- merge(imagemeta, interph[,c("name","condition","X90.","X99.","include")])
dat <- dat[dat$include=="YES",]
dat$repr_hole <- dat[,grep("X",colnames(dat))[2]] / dat[,grep("X",colnames(dat))[1]]    #dat$X90./dat$X40.
ggp <- ggplot(dat, aes(condition, repr_hole)) +
  geom_point(size=3) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "Condition", y = "pheno")
print(ggp)
t.test(
  dat$repr_hole[dat$include=="YES" & dat$condition=="+/+"],
  dat$repr_hole[dat$include=="YES" & dat$condition=="flox/flox"])
### p-val is 0.005


##Pull out example images of high and low
dat[which.min(dat$repr_hole),] #+/+
dat[which.max(dat$repr_hole),] #flox/flox


 
###################
###################
###################
###################
###################




############ Plot #holes/total area. Should be lower for flox/flox --------- meah.
dat <- merge(imagemeta, interph[,c("name","condition","nholes","include")])
dat <- dat[dat$include=="YES",]
dat$rel_holes <- dat$nholes/(dat$sum_nuc+dat$sum_cyto)
ggp <- ggplot(dat, aes(condition, rel_holes)) +
  geom_point(size=3) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "Condition", y = "pheno")
print(ggp)


############ Plot relative nuclei/cytoplasm. Might be higher in flox/flox
dat <- merge(imagemeta, interph[,c("name","condition","nholes","include")])
dat <- dat[dat$include=="YES",]
dat$rel_nuclei <- dat$sum_nuc/dat$sum_cyto
ggp <- ggplot(dat, aes(condition, rel_nuclei)) +
  geom_point(size=3) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "Condition", y = "pheno")
print(ggp)






t.test(
  interph$X80.[interph$include=="YES" & interph$color=="red"],
  interph$X80.[interph$include=="YES" & interph$color=="blue"])

sd(interph$X80.[interph$include=="YES" & interph$color=="red"])
sd(interph$X80.[interph$include=="YES" & interph$color=="blue"])


#should really average together replicates here!!!

wstat <- data.frame(
  name=interph$name,
  num=interph$num,
  stat=interph$X80)
wstat <- merge(wstat, mousecond)


wstat$pheno<-"ANY"
wstat$pheno[wstat$genotype=="flox/flox"] <- "KO"
wstat$pheno[wstat$genotype=="+/flox" & wstat$group=="a"] <- "ctrl"
wstat$pheno[wstat$genotype=="+/+"    & wstat$group=="b"] <- "ctrl"

wstat$group <- factor(wstat$group)
wstat$pheno <- factor(wstat$pheno)

themod <- lm(stat~pheno+group, wstat)
anova(themod)
