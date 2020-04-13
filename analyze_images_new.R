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
################## Read mouse annotations #############################
#######################################################################

mousecond <- read.csv("imageanalysis/allslidemeta.csv", stringsAsFactors = FALSE)
mousecond$name <- sprintf("%s-%s", mousecond$expnumber,mousecond$slidename)
#mousecond$mouse <- sprintf("imageanalysis/mydata/%s", mousecond$mouseID)


#######################################################################
################## Functions ##########################################
#######################################################################

getmask <- function(img_orig, img){
  img_diff <- img_orig[1:w,1:h,1:3]-img[1:w,1:h,1:3]
  img_diff <- img_diff[1:w,1:h,1]^2 + img_diff[1:w,1:h,2]^2 + img_diff[1:w,1:h,3]^2
  img_mask <- img_diff<0.01
  img_mask
}
# 
# 
# getcellmask <- function(img_cell){
#   #Normalize color intensities
#   n <- mean(as.double(img_orig[1:w,1:h,1]))
#   img_orig_norm <- img_orig[1:w,1:h,1]/n
#   
#   #Pull out cells
#   img_cellmask <- img_orig_norm[1:w,1:h]<0.9
#   
#   #Here could optionally fill in some gaps with the distance transformation. could then use 0.85 as a threshold
#   img_cellmask  
# }



#######################################################################
################## Classify image #####################################
#######################################################################


#outdata <- NULL
listarea <- NULL
listim <- NULL


rootdir <- "/media/mahogny/TOSHIBA EXT/rora/images"
outdir <- "/media/mahogny/TOSHIBA EXT/rora/out"
list_mice_files <- list.files(rootdir)
for(mf in list_mice_files){

  list_image_files <- list.files(sprintf("%s/%s",rootdir,mf))
  for(imf in list_image_files){

    print(sprintf("Doing %s -- %s", mf, imf))

    theim <- sprintf("%s/%s",mf,imf)
    thepath <- sprintf("%s/%s",rootdir,theim)
    
#    if(file.exists(f_orig) & file.exists(f_blood)){
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
      
    
    #### For nuclear invasion, store nuc vs cytoplasm areas  
    sum_nuc <- sum(img_class_nuc)
    sum_cyto <- sum(img_class_cyto)

    
    # kern = makeBrush(15, shape='disc')
    # img_outside <- opening(img_class_cyto, kern)
    # display(img_class_cyto)
    # display(img_outside)
    #     
    
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
    
    # #Remove the two largest areas (of which one is "0", cells, the other likely the background)
    # bg_label_area <- bg_label_area[-c(1,2)]
    # 
    # #Remove areas <3px, as these are likely technical noise
    # bg_label_area
    # sum(bg_label_area<10)
    # hist(log(bg_label_area))
    
    immeta <- data.frame(
      name=mf, 
      imname=theim, 
      sum_nuc=sum_nuc, 
      sum_cyto=sum_cyto)

    write.table(file = sprintf("%s/%s_%s.meta.csv",outdir, mf, imf), immeta, col.names = FALSE, row.names = FALSE)
    write.table(file = sprintf("%s/%s_%s.emph.csv",outdir, mf, imf), bg_label_area, col.names = FALSE, row.names = FALSE)
    
  }
}




########## Read all stats files
for(i in 1){
  listim <-rbind(listim,
                 immeta)
  listarea[[theim]] <- bg_label_area
}














allemph <- NULL
list_mice_files <- sprintf("imageanalysis/mydata/%s",list.files("imageanalysis/mydata/","mouse*"))
for(mf in list_mice_files){
  
  list_image_files <- list.files(sprintf("%s/original",mf))
  
  for(imf in list_image_files){
    
    print(sprintf("Doing %s -- %s", mf, imf))

    f_emph <- sprintf("%s/original/%s",mf,imf)
#    f_emph <- sprintf("%s/emphysema/%s",mf,imf)
    
    if(file.exists(f_emph)){
      img_emph = readImage(f_emph)

      #Normalize color intensities
      img_emph_norm <- (img_emph[1:w,1:h,1] + img_emph[1:w,1:h,2] + img_emph[1:w,1:h,3])/3
      colorMode(img_emph_norm) <- Grayscale

      #Locate background (non-cell area)
      img_outside <- img_emph_norm[1:w,1:h]>otsu(round(img_emph_norm))
      kern = makeBrush(15, shape='disc')
      img_outside <- opening(img_outside, kern)
      
      #Segment
      bg_label <- bwlabel(img_outside)
      
      #Calculate areas. Remove the two largest areas (of which one is "0", cells, the other likely the background)
      bg_label_area <- as.vector(sort(table(bg_label),decreasing = TRUE))[-c(1,2)]
      
      emph_num <- length(bg_label_area)
      emph_median <- median(bg_label_area)

      #Store summary      
      allemph <- rbind(allemph,
                       data.frame(
                         mouse=mf, 
                         image=imf, 
                         emph_len=emph_num,
                         emph_median=emph_median))
    }
  }
}



#######################################################################
################## Statistics #########################################
#######################################################################

mdata <- merge(mousecond, outdata)

ddply(mdata, .(mouseID, condition), summarise, vc=median(vesselcells))


# with distance 100
# mouseID     condition     vc
# 1  mouse_1            KO  41178
# 2 mouse_20            WT 100966
# 3 mouse_21            WT 109947
# 4 mouse_26 WT_uninfected  71046
# 5  mouse_4       control 128844


# with distance 200
# mouseID     condition       vc
# 1  mouse_1            KO 126295.5
# 2 mouse_20            WT 249788.5
# 3 mouse_21            WT 213115.5
# 4 mouse_26 WT_uninfected 153685.0
# 5  mouse_4       control 252277.0

ddply(mdata, .(mouseID, condition), summarise, vc=mean(vesselcells/vesselcircum))

# mouseID     condition        vc
# 1  mouse_1            KO 0.2405468
# 2 mouse_20            WT 0.2656210
# 3 mouse_21            WT 0.3069002
# 4 mouse_26 WT_uninfected 0.1871755
# 5  mouse_4       control 0.3153854



ddply(mdata, .(mouseID, condition), summarise, vc=mean(vesselcells/vesselarea))

# mouseID     condition         vc
# 1  mouse_1            KO 0.01416334
# 2 mouse_20            WT 0.02965755
# 3 mouse_21            WT 0.02371174
# 4 mouse_26 WT_uninfected 0.01924699
# 5  mouse_4       control 0.02117785


mdata



emphdata <- merge(mousecond, allemph)
emphdata <- emphdata[emphdata$emph_len>0,]
emphdata
ddply(emphdata, .(mouseID, condition), summarise, numareas=mean(emph_len), repr_area=mean(emph_median))

allemph



#an algorithm here!
#https://www.hindawi.com/journals/ijbi/2012/734734/
#but generally, larger areas mean emphysema 
  