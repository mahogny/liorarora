library(ggplot2)

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
######################### Load data ###################################
#######################################################################
facs_at <- read.csv("facs_adoptivetransfer/stats.csv")
#P2 all singlets... P1 AND P2?
#P9 is CD4
#P8 is Treg

facs_at$ratio_treg_cd4 <- facs_at$P8...Parent / facs_at$P9...Parent
facs_at$ratio_treg_all <- facs_at$P8...Parent 
facs_at$ratio_cd4_all  <- facs_at$P9...Parent 

facs_at_nonwt <- facs_at[facs_at$condition %in% c("rora+","rora-"),]


#######################################################################
######################### Make plots ##################################
#######################################################################

p <- ggplot(facs_at, aes(condition, ratio_treg_cd4, tissue)) + geom_point()
print(p)

p_spleen <- ggplot(facs_at_nonwt[facs_at_nonwt$tissue=="spleen",], aes(condition, ratio_treg_cd4)) + 
  geom_point(size=3) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "Spleen", y = "Treg/CD4+")
p_mln <- ggplot(facs_at_nonwt[facs_at_nonwt$tissue=="mln",], aes(condition, ratio_treg_cd4)) + 
  geom_point(size=3) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "MLN", y = "")
p_IG <- ggplot(facs_at_nonwt[facs_at_nonwt$tissue=="IG",], aes(condition, ratio_treg_cd4)) + 
  geom_point(size=3) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x = "IG", y = "")
pdf("out/facs_adoptivetransfer_treg_cd4.pdf", width = 7, height = 5)
multiplot(p_spleen, p_mln, p_IG, cols=3)
dev.off()

p <- ggplot(facs_at[facs_at$tissue=="spleen",], aes(condition, ratio_treg_all)) + geom_point()
print(p)

p <- ggplot(facs_at[facs_at$tissue=="spleen",], aes(condition, ratio_cd4_all)) + geom_point()
print(p)

#######################################################################
######################### Statistical testing #########################
#######################################################################

t.test(
  facs_at[facs_at$tissue=="IG" & facs_at$condition=="rora-",]$ratio_treg_cd4,
  facs_at[facs_at$tissue=="IG" & facs_at$condition=="rora+",]$ratio_treg_cd4)
#0.06834

t.test(
  facs_at[facs_at$tissue=="mln" & facs_at$condition=="rora-",]$ratio_treg_cd4,
  facs_at[facs_at$tissue=="mln" & facs_at$condition=="rora+",]$ratio_treg_cd4)
#0.7523

t.test(
  facs_at[facs_at$tissue=="spleen" & facs_at$condition=="rora-",]$ratio_treg_cd4,
  facs_at[facs_at$tissue=="spleen" & facs_at$condition=="rora+",]$ratio_treg_cd4)
#0.2915










