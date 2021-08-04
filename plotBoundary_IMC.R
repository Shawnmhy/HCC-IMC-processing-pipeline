library(bioimagetools)
library(jpeg)
library(tiff)
library(grid)
library(ggplot2)
library(gginnards)
library(OpenImageR)

setwd("D:/DP/Projects/HCC/Compartments")

IMC_22 <- readTIFF('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/TIFFs/Core22/Original.tiff', native = TRUE, convert = TRUE)
#bioimagetools::img(IMC_22, z=1 ,col="rgb")
IMC_8 <- readTIFF('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/TIFFs/Core8/Original.tiff', native = TRUE, convert = TRUE)


writeJPEG(IMC_22, target = "./Compartment_TIFF/Converted.jpeg", quality = 1)
writeJPEG(IMC_8, target = "./Compartment_TIFF/Converted_8.jpeg", quality = 1)

# read JPEG file

GranB_br <- readJPEG('./Compartment_TIFF/Converted.jpeg')
Original_8 <- readJPEG('./Compartment_TIFF/Converted_8.jpeg')

GranB_layer <- ggplot() +
  annotation_custom(rasterGrob(Original_8,
                               width = unit(1, "npc"), interpolate = FALSE), 
                    0, 753, 0, 759)


# reorient boundary to image




align <-  move_layers(GranB_layer, idx = 1L, position = 'bottom') +
  theme_void() +
  #geom_path(aes(Prt1.edges.conv[,1], Prt1.edges.conv[,2]), color = 'white') + , group = bdry_22[,3]
  geom_path(aes(x = bdry_8[,1], y = bdry_8[,2]), color = 'white', linetype = 'dashed') +
  xlim(0, 753) +
  ylim(0, 759) +
  coord_equal(ratio = 1)
align
ggsave(align, file= "D:/DP/Projects/HCC/Figures/aligned_TIFF_original8.png", width = 7.39, height = 7.72, units = "in", dpi = 300)







#----------- Boundary + Markers -------------#
#Core <- 22
#Marker <- 'Ki67'
#Region <- 'Immune'



#--------------- Original TIFFs to JPEG ------------------#
for(Core in c(6,7,8,9,22,23,27,28)){
  
      Original_img <- readTIFF(paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/TIFFs/Core', Core, '/Original.tiff'), native = TRUE, convert = TRUE)
      #bioimagetools::img(IMC_22, z=1 ,col="rgb")
      Original_img <- rotateFixed(Original_img, 180) %>% flipImage(mode = 'horizontal')
      
      writeJPEG(Original_img, target = paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/Converted/Core', Core, '/Original.jpeg'), quality = 1)
      
}




for(Core in c(6,7,8,9,22,23,27,28)){
  
      Core <- 8
      converted <- readJPEG(paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/Converted/Core', Core, '/Original.jpeg'))
      
      
      layer <- ggplot() +
        annotation_custom(rasterGrob(converted,
                                     width = unit(1, "npc"), interpolate = FALSE), 
                          0, 739, 0, 759)
      
      bdry <- eval(sym(paste0('bdry_', Core)))
      align <-  move_layers(layer, idx = 1L, position = 'bottom') +
        #theme_void() +
        #geom_path(aes(Prt1.edges.conv[,1], Prt1.edges.conv[,2]), color = 'white') +  group = bdry_22[,3]
        geom_path(aes(x = bdry[,1], y = bdry[,2], group = bdry$group), color = 'white', linetype = 'dashed', size = 1) +
        xlim(0, 739) +
        ylim(0, 759) +                   
        coord_equal(ratio = 1)
      
      align
      ggsave(align, file= paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/withBoundary/Core', Core, '/Original_withBoundary.jpeg'), width = 7, height = 7, units = "in", dpi = 300)
      
}



#------------------------------------------#

for(Core in c(6,7,8,9,22,23,27,28)){
  
  for(Marker in c('Ki67', 'CCR6', 'Arg1', 'PDL1', 'GranB', 'HLADR')){
    
    for(Region in c('Immune', 'Tumor')){
      Marker <- 'Lag3'
      Region <- 'Immune'
      Original_img <- readTIFF(paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/TIFFs/Core', Core, '/', Marker, '_', Region, '.tiff'), native = TRUE, convert = TRUE)
      #bioimagetools::img(IMC_22, z=1 ,col="rgb")
      Original_img <- rotateFixed(Original_img, 180) %>% flipImage(mode = 'horizontal')
      
      writeJPEG(Original_img, target = paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/Converted/Core', Core, '/', Marker, '_', Region, '.jpeg'), quality = 1)
      
    }
    
  }
}




for(Core in c(6,7,8,23,27,28)){
  
  for(Marker in c('Ki67', 'CCR6', 'Arg1', 'PDL1', 'GranB', 'HLADR')){
    
    for(Region in c('Immune')){
        
      Core <- 28
      Marker <- 'Lag3'
      Region <- 'Immune'
      
      Path <- paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/Converted/Core', Core, '/', Marker, '_', Region, '.jpeg')
      Marker_BW <- readJPEG(Path)
      
      bdry <- eval(sym(paste0('bdry_', Core)))
      Marker_Layer <- ggplot() +
        annotation_custom(rasterGrob(Marker_BW,
                                     width = unit(1, "npc"), interpolate = FALSE), 
                          0, 711, 0, 706)
      # Core 6: 749, 795
      # Core 7: 779, 780
      # Core 8: 739, 759
      # Core 9: 785, 773
      # Core 22: 739, 772
      # Core 23: 737, 746
      # Core 27: 697, 710
      # Core 28: 711, 706
      
      
      align <-  move_layers(Marker_Layer, idx = 1L, position = 'bottom') +
        theme_void() +
        #geom_point(aes(x = 0, y = 0)) +
        #geom_path(aes(Prt1.edges.conv[,1], Prt1.edges.conv[,2]), color = 'white') +  group = bdry_22[,3]
        geom_path(aes(x = bdry[,1], y = bdry[,2]), color = 'white', linetype = 'dashed', size = 2) +
        #xlim(0, 739) +
        #ylim(0, 700) +
        xlim(380, 580) +
        ylim(460, 660) + # Core 6: close: c(300,500), c(50, 250); far: c(10, 210), c(200, 400)
                         # Core 7: close: c(400,600), c(500,700); far: c(100,300), c(100,300)
                         # Core 8: close: c(250, 450), c(50, 250); far: c(200,400), c(450,650)
                         # Core 9: close: c(480, 680), c(450, 650); far: c(200,400), c(100,300)
                         # Core 23: close: c(150, 350), c(330, 530); far: c(330, 530), c(500, 700)
                         # Core 22: close: c(320, 520), c(280, 480); far: c(520, 720), c(300, 500)
                         # Core 27: close: c(200,400), c(220,420); far: c(400, 600), c(450, 650)
                         # Core 28: close: c(180, 380), c(160, 360); far: c(380, 580), c(460,660)
        coord_equal(ratio = 1)
      
      align
      ggsave(align, file= paste0('D:/DP/Projects/HCC/Figures/Compartments/Compartment_TIFF/Zoom/Core', Core, '/Far/', Marker, '_', Region, '.jpeg'), width = 7, height = 7, units = "in", dpi = 300)
      
      
    }
    
  }
}





