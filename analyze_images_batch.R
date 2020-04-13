library(EBImage)
library(CRImage)
#library(sqldf)
library(plyr)
library(stringr)

if(FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite("CRImage")
}


#######################################################################
################## Classify image #####################################
#######################################################################


rootdir <- "./images"
outdir <- "./images.out"

allmf <- list.files(rootdir) #72 of them

mf  <- allmf[Sys.getenv("mfindex")]
if(is.na(mf)){
  #For testing (if nothing given)
  mf <- "q31.16-1"
  print("Using default image")
}

#Loop over all images for sample
list_image_files <- list.files(sprintf("%s/%s",rootdir,mf))
for(imf in list_image_files){
  
  theim <- sprintf("%s/%s",mf,imf)
  thepath <- sprintf("%s/%s",rootdir,theim)
  
  if(file.exists(thepath)) {
    
    print(sprintf("Doing %s -- %s", mf, imf))
    
    img_orig = readImage(thepath)
    #    img_orig <- img_orig[1:5000,2500:5000,1:3] ###### reduce size for now
    w <- dim(img_orig)[1]
    h <- dim(img_orig)[2]
    
    #    display(img_orig[1:5000,2500:5000,1:3])
    
    
    #Classify nuc by using raw intensity
    img_class_nuc  <- img_orig[1:w,1:h,2] < 0.3
    #display(img_class_nuc)
    # 
    # #Classify nuc with light normalization
    # img_class_nuc  <- img_orig_norm[1:w,1:h,2] < 0.55
    # display(img_class_nuc)
    
    #Classify cytoplasm
    img_class_cyto <- img_orig[1:w,1:h,1] > img_orig[1:w,1:h,2] & !img_class_nuc
    display(img_class_cyto)
    display(img_class_nuc)
    
    
    #### For nuclear invasion, store nuc vs cytoplasm areas  
    sum_nuc <- sum(img_class_nuc)
    sum_cyto <- sum(img_class_cyto)
    
    
    #Join pixels together to get rid of crappy artificial holes
    kern = makeBrush(7, shape='diamond')
    img_class_cytoclean <- !opening(!img_class_cyto, kern)
    display(img_class_cytoclean)
    #display(img_class_cyto)
    #img_class_cytoclean
    
    
    #Segment
    bg_label <- bwlabel(!img_class_cytoclean)
    display(bg_label/100)
    #head(bg_label)
    
    #Calculate areas
    bg_label_area <- as.vector(sort(table(bg_label),decreasing = TRUE))
    
    immeta <- data.frame(
      name=mf, 
      imname=theim, 
      sum_nuc=sum_nuc, 
      sum_cyto=sum_cyto)
    
    write.table(file = sprintf("%s/%s_%s.meta.csv",outdir, mf, imf), immeta, col.names = FALSE, row.names = FALSE)
    write.table(file = sprintf("%s/%s_%s.emph.csv",outdir, mf, imf), bg_label_area, col.names = FALSE, row.names = FALSE)
  }
  
  
}


