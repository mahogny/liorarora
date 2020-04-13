library(EBImage)
library(CRImage)

getmask <- function(img_orig, img){
  img_diff <- img_orig[1:w,1:h,1:3]-img[1:w,1:h,1:3]
  img_diff <- img_diff[1:w,1:h,1]^2 + img_diff[1:w,1:h,2]^2 + img_diff[1:w,1:h,3]^2
  img_mask <- img_diff<0.01
  img_mask
}


#######################################################################
################## Set up pixel classifier, rgb #######################
#######################################################################

img_orig = readImage("/media/mahogny/TOSHIBA EXT/rora/classification/orig-3.png") 
#img_cell = readImage("/media/mahogny/TOSHIBA EXT/rora/classification/cytoplasm-3.png") 
img_cell = readImage("/media/mahogny/TOSHIBA EXT/rora/classification/nuc-3.png") 
img_bg   = readImage("/media/mahogny/TOSHIBA EXT/rora/classification/bg-3.png") 

w <- dim(img_orig)[1]
h<- dim(img_orig)[2]

##################################
## Normalize color intensities
norm_color <- function(img_orig){
  w <- dim(img_orig)[1]
  h<- dim(img_orig)[2]
  n <- mean(as.double(img_orig[1:w,1:h,1]))
  #n <- mean(as.double(img_orig))
  img_orig_norm <- img_orig
  for(i in 1:3){
    img_orig_norm[1:w,1:h,i] <- img_orig[1:w,1:h,i]/n
  }
  img_orig_norm
}

img_orig_norm <- norm_color(img_orig)



#Find the annotated masks
img_mask_cell <- getmask(img_orig, img_cell)
img_mask_bg   <- getmask(img_orig, img_bg)

#Extract cell and background colors
listhsv <- data.frame(
  r=as.double(img_orig_norm[1:w,1:h,1]),
  g=as.double(img_orig_norm[1:w,1:h,2]),
  b=as.double(img_orig_norm[1:w,1:h,3])
)
listhsv$rg <- listhsv$r + listhsv$g 

listhsv$div <- listhsv$r / listhsv$g 

listhsv_cell <- listhsv[as.double(img_mask_cell)==0,]
listhsv_bg   <- listhsv[as.double(img_mask_bg)==0,]

# plot(  listhsv_bg$r,   listhsv_bg$b,   cex=0.1,xlim=c(0,1),ylim=c(0,1))
# points(listhsv_cell$r, listhsv_cell$b, cex=0.1,col="red")

#### Show color distributions
plot(  listhsv_bg$r,   listhsv_bg$b,   cex=0.1,xlim=c(0.0,1.5),ylim=c(0,2))
points(listhsv_cell$r, listhsv_cell$b, cex=0.1,col="red")

plot(  listhsv_bg$g,   listhsv_bg$b,   cex=0.1,xlim=c(0.0,1.5),ylim=c(0.5,2))
points(listhsv_cell$g, listhsv_cell$b, cex=0.1,col="red")


plot(  listhsv_bg$r,   listhsv_bg$g,   cex=0.1,xlim=c(0.0,1.5),ylim=c(0.5,2))
points(listhsv_cell$r, listhsv_cell$g, cex=0.1,col="red")

plot(  listhsv_bg$rg,   listhsv_bg$b,   cex=0.1,xlim=c(0.0,3),ylim=c(0.5,2))
points(listhsv_cell$rg, listhsv_cell$b, cex=0.1,col="red")

plot(  listhsv_bg$div,   listhsv_bg$b,   cex=0.1,xlim=c(0.0,3),ylim=c(0.5,2))
points(listhsv_cell$div, listhsv_cell$b, cex=0.1,col="red")
lines(c(1,1),c(0,2))

# plot(  listhsv_bg$r,   listhsv_bg$b,   cex=0.1,xlim=c(0,1),ylim=c(0,1))
# points(listhsv_cell$r, listhsv_cell$b, cex=0.1,col="red")


##### nuclei vs bg: nuc is green<0.7

#### cytoplasm vs bg: r/g >1 ... or : r>g


#######################################################################
################## Test classifier ####################################
#######################################################################

#Classify nuc by using raw intensity
img_class_nuc  <- img_orig[1:w,1:h,2] < 0.3
display(img_class_nuc)

#Classify nuc with light normalization
img_class_nuc  <- img_orig_norm[1:w,1:h,2] < 0.55
display(img_class_nuc)

#Classify cytoplasm
img_class_cyto <- img_orig[1:w,1:h,1] > img_orig[1:w,1:h,2] & !img_class_nuc
display(img_class_cyto)
