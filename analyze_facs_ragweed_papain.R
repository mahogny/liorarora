library(ggplot2)
library(psych)


add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


readvenntable <- function(f){
  d <- read.table(f,sep=",",stringsAsFactors = FALSE)
  
  #Average over replicates
  allcond <- unique(as.character(d[1,-1]))
  m <- matrix(nrow=nrow(d)-1, ncol=length(allcond))
  for(i in 1:length(allcond)){
    v <- as.matrix(d[-1,as.character(d[1,])==as.character(allcond[i])])
    class(v) <- "double"
    m[,i] <- apply(v,1,mean)
  }
  colnames(m) <- allcond
  rownames(m) <- d[-1,1]
  list(allcond=allcond, m=m)
}



makeFacsVenn <- function(m,i){
  repRows <- function(v, n){
    m <- matrix(nrow=n,ncol=length(v))
    for(i in 1:length(v)){
      m[,i] <- v[i]
    }
    m
  }

  frac <- round(m[,i]*10)
  v <- rbind(
    repRows(c(1,0,0,0),frac[1]),
    repRows(c(0,1,0,0),frac[2]),
    repRows(c(0,0,1,0),frac[3]-frac[4]),
    repRows(c(0,0,1,1),frac[4]),
    repRows(c(0,0,0,1),frac[5]-frac[4]),
    repRows(c(0,0,0,0),frac[6])
  )
  ## Foxp3, Ifng, Il17, il13
  v  
}

makeFacsVenn_papain <- function(m,i){
  repRows <- function(v, n){
    m <- matrix(nrow=n,ncol=length(v))
    for(i in 1:length(v)){
      m[,i] <- v[i]
    }
    m
  }
  
  M<-100
  frac <- round(m[,i]*M)
  v <- rbind(
    repRows(c(1,0,0),frac[1]),
    repRows(c(0,1,0),frac[2]),
    repRows(c(0,0,0),100*M-frac[1]-frac[2])
  )
  ## il13, Foxp3
  v  
}



drawPieVenn <- function(v, 
                        ## Foxp3, Ifng, Il17, il13
                        cols=add.alpha(c("#00AA00","yellow","#FF0000","#0000CC"),0.7),
                        rads=c(1.0, 1.0, 0.95, 1.05),
                        title=""
                        ){
  plot(0,0,ylim=c(-1.5,1.5),xlim=c(-1.5,1.5),cex=0,
       xaxt="n", yaxt="n", 
       ylab="",xlab="",
       main=title)
  angs <- (1:nrow(v))*2*pi/nrow(v)

  for(curc in 1:ncol(v)){
    ni <- which.max(v[,curc]==1)+1  #there is a corner case not covered
    v[ni,curc] <- 1

    ni <- which.min(v[,curc]==1)-1  #there is a corner case not covered
    v[ni,curc] <- 1
  }
  
  allx <- cos(angs)
  ally <- sin(angs)
  polygon(allx, ally, border="black")
  for(curc in 1:ncol(v)){
    allx <- cos(angs)*v[,curc]*rads[curc]
    ally <- sin(angs)*v[,curc]*rads[curc]
    polygon(allx, ally, col=cols[curc], border = add.alpha("white",0))
  }
}



dat <- readvenntable("facs_ragweed/pie_lung.csv")
pdf("out/vennpie_ragweed_lung.pdf",w=length(dat$allcond)*3,h=3)
par(mfrow=c(1,length(dat$allcond)))
for(i in 1:length(dat$allcond)){
  drawPieVenn(makeFacsVenn(dat$m, i),title = dat$allcond[i])
}
dev.off()

dat <- readvenntable("facs_ragweed/pie_medln.csv")
pdf("out/vennpie_ragweed_medln.pdf",w=length(dat$allcond)*3,h=3)
par(mfrow=c(1,length(dat$allcond)))
for(i in 1:length(dat$allcond)){
  drawPieVenn(makeFacsVenn(dat$m, i),title = dat$allcond[i])
}
dev.off()




drawPieVenn_papain <- function(v, title){
  drawPieVenn(v, 
              ## il13, Foxp3
              cols=add.alpha(c("#0000CC", "#00AA00"),0.7),
              rads=c(1.05, 0.95),
              title = title
  )
}


dat <- readvenntable("facs_papain/pie_cd34.csv")
pdf("out/vennpie_papain_cd34.pdf",w=3*2.7,h=3*2)
par(mfrow=c(2,3))
for(i in 1:length(dat$allcond)){
  drawPieVenn_papain(makeFacsVenn_papain(dat$m, i),title = dat$allcond[i])
}
dev.off()



dat <- readvenntable("facs_papain/pie_tetraplus.csv")
pdf("out/vennpie_papain_tetraplus.pdf",w=3*2.7,h=3*2)
par(mfrow=c(2,3))
for(i in 1:length(dat$allcond)){
  drawPieVenn_papain(makeFacsVenn_papain(dat$m, i),title = dat$allcond[i])
}
dev.off()




###########################
###########################
########################### new ragweed
###########################


dat <- readvenntable("facs_ragweed/pie_lung.csv")
m <- rbind(
  data.frame(level=dat$m[1,], gene="Foxp3+", name=colnames(dat$m)),
  data.frame(level=dat$m[2,], gene="Il13+", name=colnames(dat$m)))
p<-ggplot(data=m, aes(x=name, y=level, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells (%)",fill = NULL) +
  ylim(0, 50)
p
ggsave("facs_ragweed/newout_lung.pdf")



dat <- readvenntable("facs_ragweed/pie_medln.csv")
m <- rbind(
  data.frame(level=dat$m[1,], gene="Foxp3+", name=colnames(dat$m)),
  data.frame(level=dat$m[2,], gene="Il13+", name=colnames(dat$m)))
p<-ggplot(data=m, aes(x=name, y=level, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells (%)",fill = NULL) +
  ylim(0, 50)
p
ggsave("facs_ragweed/newout_medln.pdf")




###########################
###########################
########################### new ragweed
###########################

library(ggplot2)

dat <- readvenntable("facs_papain/pie_cd34.csv")
dat
m <- rbind(
  data.frame(level=dat$m[2,], gene="Foxp3+", name=colnames(dat$m)),
  data.frame(level=dat$m[1,], gene="Il13+", name=colnames(dat$m)))
p<-ggplot(data=m, aes(x=name, y=level, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells (%)",fill = NULL) +
  ylim(0, 50)
p
ggsave("facs_papain/newout_cd34.pdf")



dat <- readvenntable("facs_papain/pie_tetraplus.csv")
m <- rbind(
  data.frame(level=dat$m[2,], gene="Foxp3+", name=colnames(dat$m)),
  data.frame(level=dat$m[1,], gene="Il13+", name=colnames(dat$m)))
p<-ggplot(data=m, aes(x=name, y=level, fill=gene)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells (%)",fill = NULL) +
  ylim(0, 50)
p
ggsave("facs_papain/newout_tetraplus.pdf")




####################################################
####################################################
####################################################


m <- data.frame(
  name=c("2W1S/Papain","Papain","PBS"),
  level=c(0.902,0.0736,0.0885)
           )

p<-ggplot(data=m, aes(x=name, y=level)) +  #, fill=gene
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "Fraction of cells (%)",fill = NULL) #+
#  ylim(0, 50)
p
#ggsave("facs_papain/newout_tetraplus.pdf")


