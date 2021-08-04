# Author: Haoyang Mi
# This script clean up the point patterns to find boundaries


library(rmarkdown)
library(dplyr)
library(reshape2)
library(ggthemes)
library(ggplot2)
library(matrixStats)
library(spatstat)
library(concaveman)
library(RANN)
library(flexclust)
library(corrplot)
library(rgeos)
library(ggpubr)
library(rstatix)
library(sp)
setwd("D:/DP/Projects/HCC/Compartments")
source('D:/DP/Projects/HCC/Functions.r')

HCCdata <- readRDS("D:/DP/Data/HCC/hccdataset")

unique(HCCdata$ctype)


# SCALE
rng <- colQuantiles(as.matrix(HCCdata[,1:32]), probs = c(0.01, 0.99))


HCCdata0 <- t((t(HCCdata[,1:32]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
HCCdata0[HCCdata0 < 0] <- 0; HCCdata0[HCCdata0 > 1] <- 1

HCCdata[,1:32] <- HCCdata0



#remove UA Noncell
HCCdata <- HCCdata[HCCdata$ctype != "UA Noncell", ]

# Include certain markers
HCCdata_ <- HCCdata[,c('CD45', 'CD16', 'CD163', 'CD68', 'CD3', 'CD4', 'CD8a', 'CD20', 'CD28', 'Foxp3', 'CD45RO', 'CD15',
                       'CCR6', 'Arg1', 'PanK', 'Ecad', 'PDL1', 'Ki67', 'HLADR', 'GranB', 'Casp3', 'Lag3', 'Collagen', 'Vimentin',
                       'Core', 'Xcoord', 'Ycoord', 'ctype', 'ctype_no')]

#cofactor <- 5
#HCCdata_[,1:24] <- asinh(HCCdata_[,1:24]/cofactor)
#HCCdata_$ctype <- HCCdata$ctype
#HCCdata_$ctype_no <- HCCdata$ctype_no

# scale using quantile
#rng <- colQuantiles(as.matrix(HCCdata_[,1:24]), probs = c(0.01, 0.99))


#HCCdata0 <- t((t(HCCdata_[,1:24]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
#HCCdata0[HCCdata0 < 0] <- 0; HCCdata0[HCCdata0 > 1] <- 1

#HCCdata_[,1:24] <- HCCdata0
#HCCdata_[,1:24] <- apply(HCCdata_[,1:24], 2, function(x){ (x - min(x)) / (max(x) - min(x)) })


#hist(as.numeric(HCCdata_$Lag3))




#------ Core 6 --------#
{
  
  Tumor1 <- fromJSON(file = './Tumor_1.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  Tumor2 <- fromJSON(file = './Tumor_2.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  Tumor3 <- fromJSON(file = './Tumor_3.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  Tumor4 <- fromJSON(file = './Tumor_4.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  Tumor5 <- fromJSON(file = './Tumor_5.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  Tumor6 <- fromJSON(file = './Tumor_6.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  ggplot() +
    geom_polygon(aes(Tumor1$X1, Tumor1$X2)) +
    geom_polygon(aes(Tumor2$X1, Tumor2$X2)) +
    geom_polygon(aes(Tumor3$X1, Tumor3$X2)) +
    geom_polygon(aes(Tumor4$X1, Tumor4$X2)) 
    #geom_polygon(aes(Tumor5$X1, Tumor5$X2)) +
    #geom_polygon(aes(Tumor6$X1, Tumor6$X2))
  
  
  # ISOLATE ISLAND: 1, 4, 6
  # 
  pts <- HCCdata_[HCCdata_$Core == '6',]
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell1 <- point.in.polygon(pts$Xcoord, pts$Ycoord, Tumor1$X1, Tumor1$X2)
  Cell2 <- point.in.polygon(pts$Xcoord, pts$Ycoord, Tumor2$X1, Tumor2$X2)
  Cell3 <- point.in.polygon(pts$Xcoord, pts$Ycoord, Tumor3$X1, Tumor3$X2)
  Cell4 <- point.in.polygon(pts$Xcoord, pts$Ycoord, Tumor4$X1, Tumor4$X2)
  Cell5 <- point.in.polygon(pts$Xcoord, pts$Ycoord, Tumor5$X1, Tumor5$X2)
  Cell6 <- point.in.polygon(pts$Xcoord, pts$Ycoord, Tumor6$X1, Tumor6$X2)
  
  inPoly1 <- which(Cell1 == 1)
  inPoly2 <- which(Cell2 == 1)
  inPoly3 <- which(Cell3 == 1)
  inPoly4 <- which(Cell4 == 1)
  inPoly5 <- which(Cell5 == 1)
  inPoly6 <- which(Cell6 == 1)
  
  # ------ HCC Compartment ---------#
  inPoly <- c(inPoly1, inPoly2, inPoly3, inPoly4, inPoly5, inPoly6)
  wrongClass1 <- which(pts$Ycoord < 55)
  wrongClass2 <- which(pts$Ycoord > 500 & pts$Xcoord < 100)
  wrongClass3 <- which(pts$Ycoord > 620 & pts$Xcoord < 300)
  wrongClass4 <- which(pts$Ycoord > 700 & pts$Xcoord < 410)
  wrongClass5 <- which(pts$Ycoord < 140 & pts$Xcoord < 100)
  
  vec <- Reduce(union, list(inPoly, wrongClass1, wrongClass2, wrongClass3, wrongClass4, wrongClass5))

  
  
  HCC_compart_6 <- pts[vec,]
  
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_6 <- pts[-vec,]
  
  
  ggplot() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    geom_path(aes(bdry_6_sec1$X1, bdry_6_sec1$X2)) +
    geom_path(aes(bdry_6_sec2$X1, bdry_6_sec2$X2)) +
    geom_path(aes(bdry_6_sec3$X1, bdry_6_sec3$X2)) +
    geom_path(aes(bdry_6_sec4$X1, bdry_6_sec4$X2)) +
    geom_path(aes(bdry_6_sec5$X1, bdry_6_sec5$X2)) +
    geom_path(aes(bdry_6_sec6$X1, bdry_6_sec6$X2)) +
    geom_point(aes(Imm_compart_6$Xcoord, Imm_compart_6$Ycoord), color = 'blue') +
    geom_point(aes(HCC_compart_6$Xcoord, HCC_compart_6$Ycoord), color = 'red') 
  
  
  
  # ----- Compute boundary -------#
  bdry_6_sec1 <- Tumor1
  bdry_6_sec2 <- Tumor4
  bdry_6_sec3 <- Tumor5
  bdry_6_sec4 <- Tumor6

  
  
  # ----- Compute boundary -------#
  bdry_6_sec5 <- Tumor2[c(2500:4363),]
  bdry_6_sec6 <- Tumor3[c(500:1600),]

  plot(pts$Xcoord, pts$Ycoord)
  points(bdry_6_sec6)
  points(bdry_6_sec5)
  points(bdry_6_sec4)
  points(bdry_6_sec1)
  points(bdry_6_sec2)
  points(bdry_6_sec3)
  
  
  
  bdry_6_sec1$group <- 1
  bdry_6_sec2$group <- 2
  bdry_6_sec3$group <- 3
  bdry_6_sec4$group <- 4
  bdry_6_sec5$group <- 5
  bdry_6_sec6$group <- 6
  
  bdry_6 <- rbind(bdry_6_sec1, bdry_6_sec2, bdry_6_sec3, bdry_6_sec4, bdry_6_sec5, bdry_6_sec6)
  
  
  bdry_6_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_6_sec1[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_6_sec3[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_6_sec2[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_6_sec4[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_6_sec5[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_6_sec6[,1:2])),1))))/1000 
  
  # compute distance to boundary
  d2bdry_Imm_6 <- dist2(Imm_compart_6[, c('Xcoord', 'Ycoord')], bdry_6[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_6)
  colnames(d2bdry_Imm_6)[1] <- 'dist'
  d2bdry_HCC_6 <- dist2(HCC_compart_6[, c('Xcoord', 'Ycoord')], bdry_6[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_6)
  colnames(d2bdry_HCC_6)[1] <- 'dist'
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_6$Xcoord/1000, d2bdry_Imm_6$Ycoord/1000, fill = d2bdry_Imm_6$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_6$Xcoord/1000, d2bdry_HCC_6$Ycoord/1000, fill = d2bdry_HCC_6$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_6_sec1$X1/1000, bdry_6_sec1$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_6_sec2$X1/1000, bdry_6_sec2$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_6_sec3$X1/1000, bdry_6_sec3$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_6_sec4$X1/1000, bdry_6_sec4$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_6_sec5$X1/1000, bdry_6_sec5$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_6_sec6$X1/1000, bdry_6_sec6$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/6_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300) 
  
}




#------ Core 7 --------#
{
  
  json_data <- fromJSON(file = './Immune_7.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  ggplot() +
    geom_polygon(aes(json_data$X1, json_data$X2))
  
  
  
  # 
  pts <- HCCdata_[HCCdata_$Core == '7',]
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, json_data$X1, json_data$X2)
  
  inPoly <- which(Cell == 1)
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_7 <- pts[inPoly,]
  
  
  # ------ HCC Compartment ---------#
  
  HCC_compart_7 <- pts[-inPoly,]
  
  
  
  # ----- Compute boundary -------#
  bdry_7 <- json_data
    
  
  bdry_7_sec1 <- bdry_7[1:1000,]
  bdry_7_sec2 <- bdry_7[1500:2361,]
  
  bdry_7 <- rbind(bdry_7_sec2, bdry_7_sec1)
  
  ggplot() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    geom_path(aes(bdry_7$X1, bdry_7$X2)) +

    geom_point(aes(Imm_compart_7$Xcoord, Imm_compart_7$Ycoord), color = 'blue') +
    geom_point(aes(HCC_compart_7$Xcoord, HCC_compart_7$Ycoord), color = 'red') 
  
  bdry_7_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_7[,1:2])),1))))/1000 
  
  
  # compute distance to boundary
  d2bdry_Imm_7 <- dist2(Imm_compart_7[, c('Xcoord', 'Ycoord')], bdry_7[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_7)
  colnames(d2bdry_Imm_7)[1] <- 'dist'
  d2bdry_HCC_7 <- dist2(HCC_compart_7[, c('Xcoord', 'Ycoord')], bdry_7[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_7)
  colnames(d2bdry_HCC_7)[1] <- 'dist'
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_7$Xcoord/1000, d2bdry_Imm_7$Ycoord/1000, fill = d2bdry_Imm_7$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_7$Xcoord/1000, d2bdry_HCC_7$Ycoord/1000, fill = d2bdry_HCC_7$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_7_sec1$X1/1000, bdry_7_sec1$X2/1000), color = '#03f67d', size = 2) +
    geom_line(aes(bdry_7_sec2$X1/1000, bdry_7_sec2$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/7_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
    
  
}


#------ Core 8 --------#
{
  
  json_data <- fromJSON(file = './Tumor_8.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  ggplot() +
    geom_polygon(aes(json_data$X1, json_data$X2))
  
  
  
  # 
  pts <- HCCdata_[HCCdata_$Core == '8',]
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, json_data$X1, json_data$X2)
  
  inPoly <- which(Cell == 1)
  wrongClass <- which(pts$Ycoord > 500)
  
  vec <- union(inPoly, wrongClass)
  # ------ HCC Compartment ---------#
  
  HCC_compart_8 <- pts[vec,]
  
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_8 <- pts[-vec,] 
  

  
  
  # ----- Compute boundary -------#
  bdry_8_sec1 <- json_data[c(3000:4574),]
  bdry_8_sec2 <- json_data[c(1:200),]
  
  plot(pts$Xcoord, pts$Ycoord)
  #points(bdry_8)
  points(bdry_8_sec2)
  
  
  bdry_8_sec1$group <- 1
  bdry_8_sec2$group <- 2
  
  bdry_8 <- rbind(bdry_8_sec1, bdry_8_sec2)
  
  bdry_8_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_8_sec1[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_8_sec2[,1:2])),1))))/1000 
  
  
  # compute distance to boundary
  d2bdry_Imm_8 <- dist2(Imm_compart_8[, c('Xcoord', 'Ycoord')], bdry_8[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_8)
  colnames(d2bdry_Imm_8)[1] <- 'dist'
  d2bdry_HCC_8 <- dist2(HCC_compart_8[, c('Xcoord', 'Ycoord')], bdry_8[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_8)
  colnames(d2bdry_HCC_8)[1] <- 'dist'
  
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_8$Xcoord/1000, d2bdry_Imm_8$Ycoord/1000, fill = d2bdry_Imm_8$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_8$Xcoord/1000, d2bdry_HCC_8$Ycoord/1000, fill = d2bdry_HCC_8$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_8_sec1$X1/1000, bdry_8_sec1$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_8_sec2$X1/1000, bdry_8_sec2$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/8_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
  
  
  
}


#------ Core 9 --------#
{
  
  json_data <- fromJSON(file = './Immune_9.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  ggplot() +
    geom_polygon(aes(json_data$X1, json_data$X2))
  
  
  
  # 
  pts <- HCCdata_[HCCdata_$Core == '9',]
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, json_data$X1, json_data$X2)
  
  inPoly <- which(Cell == 1)
  wrongClass <- which(pts$Xcoord > 600 & pts$Ycoord > 520)
  
  vec <- union(inPoly, wrongClass)
  # ------ HCC Compartment ---------#
  
  HCC_compart_9 <- pts[-vec,]
  
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_9 <- pts[vec,] 
  
  
  
  
  # ----- Compute boundary -------#
  bdry_9_sec1 <- json_data[c(1500:2339),]
  bdry_9_sec2 <- json_data[c(1:950),]
  
  plot(pts$Xcoord, pts$Ycoord)
  #points(bdry_8)
  points(bdry_9_sec1)
  points(bdry_9_sec2)
  
  
  bdry_9_sec1$group <- 1
  bdry_9_sec2$group <- 2
  
  bdry_9 <- rbind(bdry_9_sec1, bdry_9_sec2)
  
  bdry_9_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_9_sec1[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_9_sec2[,1:2])),1))))/1000 
  
  
  # compute distance to boundary
  d2bdry_Imm_9 <- dist2(Imm_compart_9[, c('Xcoord', 'Ycoord')], bdry_9[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_9)
  colnames(d2bdry_Imm_9)[1] <- 'dist'
  d2bdry_HCC_9 <- dist2(HCC_compart_9[, c('Xcoord', 'Ycoord')], bdry_9[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_9)
  colnames(d2bdry_HCC_9)[1] <- 'dist'
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_9$Xcoord/1000, d2bdry_Imm_9$Ycoord/1000, fill = d2bdry_Imm_9$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_9$Xcoord/1000, d2bdry_HCC_9$Ycoord/1000, fill = d2bdry_HCC_9$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_9_sec1$X1/1000, bdry_9_sec1$X2/1000), color = '#03f67d', size = 2) +
    geom_path(aes(bdry_9_sec2$X1/1000, bdry_9_sec2$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/9_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
  
  
}



#------ Core 22 --------#
{
  
  json_data <- fromJSON(file = './Tumor_22.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  
  
  
  # 
  pts <- HCCdata_[HCCdata_$Core == '22',]
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, json_data$X1, json_data$X2)
  
  inPoly <- which(Cell == 1)
  wrongClass <- which(pts$Ycoord < 200 & pts$Xcoord > 410)
  wrongClass1 <- which(pts$Xcoord > 610)
  wrongClass2 <- which(pts$Xcoord < 300 & pts$Ycoord  > 700)
  wrongClass3 <- which(pts$Xcoord < 30 & pts$Ycoord > 360)
  
  vec <- Reduce(union, list(inPoly, wrongClass, wrongClass1, wrongClass2, wrongClass3))
  
  # ------ HCC Compartment ---------#
  
  HCC_compart_22 <- pts[vec,]
  
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_22 <- pts[-vec,]
  
  
  # ----- Compute boundary -------#
  bdry_22_sec1 <- json_data[c(1:2200),]
  plot(pts$Xcoord, pts$Ycoord)
  points(bdry_22_sec1)
  
  
  bdry_22_sec2 <- json_data[c(2800: 3300),]
  
  plot(pts$Xcoord, pts$Ycoord)
  points(bdry_22_sec2)
  
  
  
  
  
  bdry_22_sec1$group <- 1
  bdry_22_sec2$group <- 2
  bdry_22 <- rbind(bdry_22_sec1, bdry_22_sec2)
  
  bdry_22_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_22_sec1[,1:2])),1))))/1000 +
    gLength(SpatialLines(list(Lines(list(Line(bdry_22_sec2[,1:2])),1))))/1000 
  
  
  # compute distance to boundary
  d2bdry_Imm_22 <- dist2(Imm_compart_22[, c('Xcoord', 'Ycoord')], bdry_22[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_22)
  colnames(d2bdry_Imm_22)[1] <- 'dist'
  d2bdry_HCC_22 <- dist2(HCC_compart_22[, c('Xcoord', 'Ycoord')], bdry_22[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_22)
  colnames(d2bdry_HCC_22)[1] <- 'dist'
  
  ggplot() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    geom_path(aes(bdry_22_sec1$X1, bdry_22_sec1$X2)) +
    geom_point(aes(Imm_compart_22$Xcoord, Imm_compart_22$Ycoord), color = 'BLUE') +
    geom_point(aes(HCC_compart_22$Xcoord, HCC_compart_22$Ycoord), color = 'red') 
  
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_22$Xcoord/1000, d2bdry_Imm_22$Ycoord/1000, fill = d2bdry_Imm_22$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_22$Xcoord/1000, d2bdry_HCC_22$Ycoord/1000, fill = d2bdry_HCC_22$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_22_sec1$X1/1000, bdry_22_sec1$X2/1000), color = '#03f67d', size = 2) +
    geom_line(aes(bdry_22_sec2$X1/1000, bdry_22_sec2$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
    #scale_fill_gradientn(colours = pal) +
  
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/22_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
  
  
}


#------ Core 23 --------#
{
  
  json_data <- fromJSON(file = './Tumor_23.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  plot(json_data)
  
  
  # 
  pts <- HCCdata_[HCCdata_$Core == '23',]
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, json_data$X1, json_data$X2)
  
  inPoly <- which(Cell == 1)
  wrongClass <- which(pts$Ycoord > 700)
  wrongClass1 <- which(pts$Xcoord > 600 & pts$Ycoord < 300)
  wrongClass2 <- which(pts$Ycoord < 40)
  wrongClass3 <- which(pts$Xcoord < 60 & pts$Ycoord > 500)
  wrongClass3 <- which(pts$Xcoord > 400)
  
  vec <- Reduce(union, list(inPoly, wrongClass, wrongClass1, wrongClass2, wrongClass3))
  # ------ HCC Compartment ---------#
  
  HCC_compart_23 <- pts[vec,]
  
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_23 <- pts[-vec,]
  
  
  
  ggplot() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    geom_polygon(aes(json_data$X1, json_data$X2)) +
    geom_point(aes(Imm_compart_23$Xcoord, Imm_compart_23$Ycoord), color = 'blue') +
    geom_point(aes(HCC_compart_23$Xcoord, HCC_compart_23$Ycoord), color = 'red')
  
  
  
  
  
  # ----- Compute boundary -------#
  bdry_23_sec1 <- json_data[c(280:2500),]
  plot(pts$Xcoord, pts$Ycoord)
  points(bdry_23_sec1)
  
  
  #bdry_23_sec2 <- json_data[c(3900:4100),]
  
  #plot(pts$Xcoord, pts$Ycoord)
  #points(bdry_23_sec2)
  
  
  
  #bdry_23_sec1$group <- 1

  bdry_23 <- bdry_23_sec1
  
  bdry_23_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_23[,1:2])),1))))/1000
  
  # compute distance to boundary
  d2bdry_Imm_23 <- dist2(Imm_compart_23[, c('Xcoord', 'Ycoord')], bdry_23[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_23)
  colnames(d2bdry_Imm_23)[1] <- 'dist'
  d2bdry_HCC_23 <- dist2(HCC_compart_23[, c('Xcoord', 'Ycoord')], bdry_23[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_23)
  colnames(d2bdry_HCC_23)[1] <- 'dist'
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_23$Xcoord/1000, d2bdry_Imm_23$Ycoord/1000, fill = d2bdry_Imm_23$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_23$Xcoord/1000, d2bdry_HCC_23$Ycoord/1000, fill = d2bdry_HCC_23$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_23_sec1$X1/1000, bdry_23_sec1$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/23_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
  
}

#------ Core 27 --------#
{
  
  json_data <- fromJSON(file = './Immune_27.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  plot(json_data)
  
  
  # 
  pts <- HCCdata_[HCCdata_$Core == '27',]
  
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, json_data$X1, json_data$X2)
  
  inPoly <- which(Cell == 1)
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_27 <- pts[inPoly,]
  
  
  # ------ HCC Compartment ---------#
  
  HCC_compart_27 <- pts[-inPoly,]
  
  
  
  # ----- Compute boundary -------#
  bdry_27 <- json_data
  
  
  bdry_27_L <- 
    gLength(SpatialLines(list(Lines(list(Line(json_data[,1:2])),1))))/1000 
  
  # compute distance to boundary
  d2bdry_Imm_27 <- dist2(Imm_compart_27[, c('Xcoord', 'Ycoord')], bdry_27[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_27)
  colnames(d2bdry_Imm_27)[1] <- 'dist'
  d2bdry_HCC_27 <- dist2(HCC_compart_27[, c('Xcoord', 'Ycoord')], bdry_27[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_27)
  colnames(d2bdry_HCC_27)[1] <- 'dist'
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_27$Xcoord/1000, d2bdry_Imm_27$Ycoord/1000, fill = d2bdry_Imm_27$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_27$Xcoord/1000, d2bdry_HCC_27$Ycoord/1000, fill = d2bdry_HCC_27$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_27$X1/1000, bdry_27$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/27_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
  
}

#------ Core 28 --------#
{
  
  json_data <- fromJSON(file = './Immune_28.json')$geometry$coordinates %>% data.frame() %>% t() %>% data.frame()
  #plot(json_data)
  
  
  ggplot() +
    geom_polygon(aes(json_data$X1, json_data$X2))
  # 
  pts <- HCCdata_[HCCdata_$Core == '28',]
  
  
  
  
  # ----- Compute boundary -------#
  bdry_28 <- json_data[c(4420:5180),]
  plot(pts$Xcoord, pts$Ycoord)
  points(bdry_28)
  
  #ggplot() +
  #  geom_polygon(aes(bdry_28$X1, bdry_28$X2))
  
  
  # ----- Cells in Tumor Compartment ------#
  
  Cell <- point.in.polygon(pts$Xcoord, pts$Ycoord, bdry_28$X1, bdry_28$X2)
  
  inPoly <- which(Cell == 1)
  
  # ------ HCC Compartment ---------#
  
  Imm_compart_28 <- pts[inPoly,]
  
  
  # ------ HCC Compartment ---------#
  
  HCC_compart_28 <- pts[-inPoly,]
  
  
  
  # ----- Compute boundary -------#
  
  bdry_28_L <- 
    gLength(SpatialLines(list(Lines(list(Line(bdry_28[,1:2])),1))))/1000 
  
  # compute distance to boundary
  d2bdry_Imm_28 <- dist2(Imm_compart_28[, c('Xcoord', 'Ycoord')], bdry_28[,1:2]) %>% apply(1, min) %>% cbind(Imm_compart_28)
  colnames(d2bdry_Imm_28)[1] <- 'dist'
  d2bdry_HCC_28 <- dist2(HCC_compart_28[, c('Xcoord', 'Ycoord')], bdry_28[,1:2]) %>% apply(1, min) %>% cbind(HCC_compart_28)
  colnames(d2bdry_HCC_28)[1] <- 'dist'
  
  p <- ggplot() +
    theme_bw() +
    #geom_point(aes(pts$Xcoord, pts$Ycoord), color = 'red') +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 32),
          legend.position = 'none') +
    geom_point(aes(d2bdry_Imm_28$Xcoord/1000, d2bdry_Imm_28$Ycoord/1000, fill = d2bdry_Imm_28$dist), size = 4, shape = 21) +
    geom_point(aes(d2bdry_HCC_28$Xcoord/1000, d2bdry_HCC_28$Ycoord/1000, fill = d2bdry_HCC_28$dist), size = 4, shape = 21) +
    geom_path(aes(bdry_28$X1/1000, bdry_28$X2/1000), color = '#03f67d', size = 2) +
    xlab('x, mm') +
    ylab('y, mm') +
    scale_fill_gradient(low = "#82bed7", high = "#fe0101")
  #scale_fill_gradientn(colours = pal) +
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/28_comparment_border.png"), width = 8, height = 8, units = "in", dpi = 300)
  
}


Imm_in_Tumor_all <- data.frame(matrix(nrow = 0, ncol = 0))
Tumor_in_Imm_all <- data.frame(matrix(nrow = 0, ncol = 0))
simple_Imm_in_Tumor_all <- data.frame(matrix(nrow = 0, ncol = 0))
simple_Tumor_in_Imm_all <- data.frame(matrix(nrow = 0, ncol = 0))
profile_Ki67_all_cores <- data.frame(matrix(nrow = 0, ncol = 0))

# read compartment area data
compartArea <- read.csv('compartment_areas.csv')

for(core in c(6, 7, 8, 9, 22, 23, 27, 28)){
  
  # core data
  #core <- 28
  
  coreData <- HCCdata[HCCdata$Core == core, ]
  
  response <- unique(coreData$response)
  # boundary data
  bdry <- eval(sym(paste0('bdry_', core)))
  
  
  
  # get IMM compartment cells
  Imm_compart <- eval(sym(paste0('d2bdry_Imm_', core)))
  
  # get HCC compartment cells
  HCC_compart <- eval(sym(paste0('d2bdry_HCC_', core)))
  
  
  # maximum distance for Imm compatment
  maxImm <- max(Imm_compart$dist)
  
  # HCC distance for Imm compatment
  maxHCC <- max(HCC_compart$dist)
  
  # get cell types with non-zero counts
  
  cellAval <- as.character(data.frame(table(coreData$ctype_no))$Var1)
  
  
  # boundary length of the corresponding core
  bdry_length <- eval(sym(paste0('bdry_', core, '_L')))
  
  
  
  # profile data frame
  profile <- data.frame(matrix(nrow = 0, ncol = 0))
  profile1 <- data.frame(matrix(nrow = 0, ncol = 0))
  
  profile_Ki67 <- data.frame(matrix(nrow = 0, ncol = 0))
  profile1_Ki67 <- data.frame(matrix(nrow = 0, ncol = 0))
  
  for(ctype in unique(cellAval)){
    
    
    cell_cype <- Imm_compart[Imm_compart$ctype_no == ctype,]
    
    for(band in seq(0, ceiling(maxImm), 20)/2){
      
      
      count <- nrow(cell_cype[cell_cype$dist >= band & cell_cype$dist < band + 20,]) / bdry_length
      

      profile <- rbind(profile, cbind(band, count, ctype))
      
      
    }
    
    cell_cype <- HCC_compart[HCC_compart$ctype_no == ctype,]
    
    for(band in seq(0, ceiling(maxHCC), 20)/2){
      
      
      count <- nrow(cell_cype[cell_cype$dist >= band & cell_cype$dist < band + 20,]) / bdry_length
      
      
      profile1 <- rbind(profile1, cbind(band, count, ctype))
      
      
    } 
    
  }
  
  
  #----- Ki67 gradient ----#
  stepsize <- 20
  for(band in seq(0, ceiling(maxImm), stepsize)/2){
    
    Ki67_mean <- Imm_compart[Imm_compart$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),] %>%
      filter(dist >= band & dist < band + stepsize) %>%
      summarise(Mean = mean(Ki67)) %>%
      as.numeric()
    
    
    profile_Ki67 <- rbind(profile_Ki67, cbind(band, Ki67_mean))
    
    
  }
  
  for(band in seq(0, ceiling(maxHCC), stepsize)/2){
    
    Ki67_mean <- HCC_compart[HCC_compart$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),] %>%
      filter(dist >= band & dist < band + stepsize) %>%
      summarise(Mean = mean(Ki67)) %>%
      as.numeric()
    
    profile1_Ki67 <- rbind(profile1_Ki67, cbind(band, Ki67_mean))
    
    
  } 
  
  #----------#
  
  
  profile$band <- as.numeric(as.character(profile$band)) + 10
  profile$count <- as.numeric(as.character(profile$count))
  
  profile1$band <- 0 - (as.numeric(as.character(profile1$band)) + 10)
  profile1$count <- as.numeric(as.character(profile1$count))
  profile_all <- rbind(profile, profile1)
  
  
  
  profile_Ki67$band <- as.numeric(as.character(profile_Ki67$band)) + 10
  profile_Ki67$Ki67_mean <- as.numeric(as.character(profile_Ki67$Ki67_mean))
  
  profile1_Ki67$band <- 0 - (as.numeric(as.character(profile1_Ki67$band)) + 10)
  profile1_Ki67$Ki67_mean <- as.numeric(as.character(profile1_Ki67$Ki67_mean))
  
  profile_Ki67_all <- rbind(profile_Ki67, profile1_Ki67) %>%
    mutate(core = core)
  
  profile_Ki67_all_cores <- rbind(profile_Ki67_all_cores, profile_Ki67_all)
  # calculate distance of coredate to boundary
  p <- ggplot(profile_all, aes(band/1000, count, group = ctype, color = as.character(ctype), shape = ctype)) +
    theme_bw() +
    geom_line(size = 1) +
    geom_point(size = 3) +
    #scale_fill_gradientn(colours = pal) +
    theme(legend.position = 'none',
          axis.title = element_text(size = 24), 
          axis.text = element_text(size = 22)) +
    geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
    scale_shape_manual(values = c('ctype1' = 1, 'ctype2' = 2, 'ctype3' = 3,
                                  'ctype5' = 4, 'ctype6' = 5, 'ctype7' = 6,
                                  'ctype8' = 7, 'ctype9' = 8, 'ctype10' = 9,
                                  'ctype11' = 10, 'ctype12' = 11, 'ctype13' = 12,
                                  'ctype14' = 13, 'ctype15' = 14, 'ctype16' = 15,
                                  'ctype17' = 16, 'ctype18' = 17)) +
    scale_color_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                  'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                  'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                  'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                  'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                  'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
    xlab('Distance to boundary, mm') +
    ylab('Normalized Counts (/mm)') 
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Profiles/", core,"_normalized.png"), width = 25, height = 6, units = "in", dpi = 300)
  
  
  
  # glm
  
  #profile_all_wide <- reshape(profile_all, idvar = "band", timevar = "ctype", direction = "wide")
  
  

  #wideTable <- dcast(data = profile_all, formula = band~ctype, value.var = 'count')
  
  #cnname <- as.character(profile_all[profile_all$ctype %in% HCCdata$ctype,]$ctype)
  #wideTable <- wideTable[,-6]
  #M <- cor(wideTable[, 2:13])
  
  
  # matrix of the p-value of the correlation
  #p.mat <- cor.mtest(wideTable[, 2:13])
  
  
  
  #col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  #corrplot(M, method="color", col=col(200),  
  #         type="upper", order="hclust", 
  #         addCoef.col = "black", # Add coefficient of correlation
  #         tl.col="black", tl.srt=45, #Text label color and rotation
           # Combine with significance
  #         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
  #         # hide correlation coefficient on the principal diagonal
  #         diag=FALSE 
  #)
  
  # Immune compartment area for the current core
  ImmArea <- compartArea[compartArea$core == core, 'areaImm.um'] %>%
    as.character() %>%
    as.numeric()
  HCCArea <- compartArea[compartArea$core == core, 'areaHCC.um'] %>%
    as.character() %>%
    as.numeric()
  
  # Immune cell in tumor
  Imm_in_Tumor <- HCC_compart %>%
    filter(ctype_no %in% c('ctype1', 'ctype2', 'ctype3', 'ctype5', 'ctype7', 'ctype8', 'ctype10', 'ctype11', 'ctype13', 'ctype15', 'ctype18')) %>%
    with(table(ctype_no)) %>%
    data.frame() %>%
    mutate(density = Freq / HCCArea)
  
  # Tumor cell in immune
  Tumor_in_Imm <- Imm_compart %>%
    filter(ctype_no %in% c('ctype9', 'ctype12', 'ctype14', 'ctype16', 'ctype17')) %>%
    with(table(ctype_no)) %>%
    data.frame() %>%
    mutate(density = Freq / ImmArea)
  
  
  
  Imm_in_Tumor_all <- rbind(Imm_in_Tumor_all, cbind(core, Imm_in_Tumor, response))
  
  Tumor_in_Imm_all <- rbind(Tumor_in_Imm_all, cbind(core, Tumor_in_Imm, response))
  
  
  Imm_in_Tumor <- HCC_compart %>%
    filter(ctype_no %in% c('ctype1', 'ctype2', 'ctype3', 'ctype5', 'ctype7', 'ctype8', 'ctype10', 'ctype11', 'ctype13', 'ctype15', 'ctype18')) %>%
    nrow()
  
  # Tumor cell in immune
  Tumor_in_Imm <- Imm_compart %>%
    filter(ctype_no %in% c('ctype9', 'ctype12', 'ctype14', 'ctype16', 'ctype17')) %>%
    nrow()

  
  simple_Imm_in_Tumor_all <- rbind(simple_Imm_in_Tumor_all, cbind(core, Imm_in_Tumor / HCCArea, response))
  
  simple_Tumor_in_Imm_all <- rbind(simple_Tumor_in_Imm_all, cbind(core, Tumor_in_Imm / ImmArea, response))
  
}


profile_Ki67_all_cores$core <- as.character(profile_Ki67_all_cores$core)


####
####
#### group to R and NR
{
  R <- profile_Ki67_all_cores %>%
    filter(core %in% c(22, 23))
  
  NR <- profile_Ki67_all_cores %>%
    filter(core %in% c(6, 7, 8, 9, 27, 28))
  
  reCalc <- data.frame(matrix(nrow = 0, ncol = 0))
  
  for(flag in c('R', 'NR')){
    
    dat <- eval(sym(flag))
    
    
    for(dist in unique(dat$band)){
      reKi67 <- dat %>%
        filter(band == dist) %>%
        summarise(Mean = mean(Ki67_mean)) %>%
        as.numeric()    
      
      reCalc <- rbind(reCalc, cbind(dist, reKi67, flag))
    }
  }
  
  
  reCalc$reKi67 <- as.numeric(as.character(reCalc$reKi67))
  reCalc$dist <- as.numeric(as.character(reCalc$dist))
  
  p <- ggplot(reCalc, aes(dist/1000, reKi67, group = flag, color = as.character(flag))) +
    theme_bw() +
    geom_line(size = 1) +
    geom_point(size = 3) +
    #scale_fill_gradientn(colours = pal) +
    theme(legend.position = 'none',
          axis.title = element_text(size = 24), 
          axis.text = element_text(size = 22)) +
    geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
    xlab('Distance to boundary, mm') +
    ylab('Ki67 Intensity') 
  p
  
  
}






simple_Imm_in_Tumor_all$V2 <- as.numeric(as.character(simple_Imm_in_Tumor_all$V2))
colnames(simple_Imm_in_Tumor_all)[2] <- 'density'
simple_Tumor_in_Imm_all$V2 <- as.numeric(as.character(simple_Tumor_in_Imm_all$V2))
colnames(simple_Tumor_in_Imm_all)[2] <- 'density' 
# marker distribution matrix

NR_data <- rbind(d2bdry_HCC_6, d2bdry_HCC_7, d2bdry_HCC_8, d2bdry_HCC_9, d2bdry_HCC_27, d2bdry_HCC_28,  d2bdry_Imm_7, d2bdry_Imm_8,
                 d2bdry_Imm_9, d2bdry_Imm_6, d2bdry_Imm_27, d2bdry_Imm_28)

R_data <- rbind(d2bdry_HCC_22, d2bdry_HCC_23, d2bdry_Imm_22, d2bdry_Imm_23)

markerGradient_R <- data.frame(matrix(nrow = 0, ncol = 3))
markerGradient_NR <- data.frame(matrix(nrow = 0, ncol = 3))

for(col in c('CCR6', 'Arg1', 'PDL1', 'HLADR', 'GranB', 'Ki67')){
  
  #col <- 'PDL1'
  R_data <- rbind(d2bdry_HCC_22, d2bdry_HCC_23, d2bdry_Imm_22, d2bdry_Imm_23)
  NR_data <- rbind(d2bdry_HCC_6, d2bdry_HCC_7, d2bdry_HCC_8, d2bdry_HCC_9, d2bdry_HCC_27, d2bdry_HCC_28,  d2bdry_Imm_7, d2bdry_Imm_8,
                   d2bdry_Imm_9, d2bdry_Imm_6, d2bdry_Imm_27, d2bdry_Imm_28)
  
  
  
  if(col == 'GranB' | col == 'Lag3'){
    R_data <- R_data[R_data$ctype_no == 'ctype8',]
    NR_data <- NR_data[NR_data$ctype_no == 'ctype8',]
    
  }
  
  if(col == 'Ki67' | col == 'PDL1' | col == 'CCR6'){
    R_data <- R_data[R_data$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
    NR_data <- NR_data[NR_data$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
    
  }
  
  if(col == 'Arg1'){
    R_data <- R_data[R_data$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),]
    NR_data <- NR_data[NR_data$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),]
    
  }
  
  
  
  close <- data.frame(R_data[R_data$dist < 40 , col]) %>%
    mutate(class = 'close', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  far <- data.frame(R_data[R_data$dist >= 40, col]) %>%
    mutate(class = 'far', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  markerGradient_R <- rbind(markerGradient_R, rbind(close, far))
  
  
  temp <- data.frame(rbind(close, far))
  
  close <- data.frame(NR_data[NR_data$dist < 40 , col]) %>%
    mutate(class = 'close', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  far <- data.frame(NR_data[NR_data$dist >= 40 , col]) %>%
    mutate(class = 'far', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  markerGradient_NR <- rbind(markerGradient_NR, rbind(close, far))
  
  
  
}




markerGradient_R$response <- 'R'
markerGradient_NR$response <- 'NR'

markerGradient <- rbind(markerGradient_R, markerGradient_NR)


stat.test <- markerGradient %>%
  group_by(response, marker) %>%
  wilcox_test(value~class) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()



stat.test <- stat.test %>% add_xy_position(x = "class")


jpeg('./Figures/Functional_Marker_Profile_Outcome.jpeg', units="in", width=15, height=8, res=300)

# adjusted p value

ggboxplot(markerGradient, x = 'class', y = 'value', palette = 'jco', color = 'class', facet.by = c('response', 'marker')) +
  #stat_boxplot(geom='errorbar', linetype=1, width=0.5)+  #whiskers
  #geom_boxplot(outlier.shape=1)+    
  #theme_bw() +
  #stat_compare_means(aes(label=..p.adj..), comparisons = list(c('close', 'far')), method = 'wilcox.test', label.y = max(markerGradient_R$value) + 0.05,
  #                   p.adjust.methods = 'fdr', size = 5) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = max(markerGradient$value) + 0.05, label.size = 6) +

  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 16)) +
  ylab('Normalized Expression') +
  ylim(0, 1.1)
#scale_y_continuous(breaks = seq(0, 1.2, by = 0.25)) 

#facet_grid(vars(marker), var(response))
dev.off()



stat.test_prox <- markerGradient %>%
  group_by(class, marker) %>%
  wilcox_test(value~response) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()



stat.test_prox <- stat.test_prox %>% add_xy_position(x = "response")


jpeg('./Figures/Functional_Marker_Profile_Prox.jpeg', units="in", width=15, height=8, res=300)

ggboxplot(markerGradient, x = 'response', y = 'value', palette = 'lancet', color = 'response', facet.by = c('class', 'marker')) +
  #stat_boxplot(geom='errorbar', linetype=1, width=0.5)+  #whiskers
  #geom_boxplot(outlier.shape=1)+    
  #theme_bw() +
  #stat_compare_means(aes(label=..p.adj..), comparisons = list(c('close', 'far')), method = 'wilcox.test', label.y = max(markerGradient_R$value) + 0.05,
  #                   p.adjust.methods = 'fdr', size = 5) +
  stat_pvalue_manual(stat.test_prox, label = "p.adj.signif", y.position = max(markerGradient$value) + 0.05, label.size = 6) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 16)) +
  ylab('Normalized Expression') +
  ylim(0, 1.1)
dev.off()






# ----------------------------------------------- #
#       Split to Immune and HCC compartment       #
# ----------------------------------------------- #






# marker distribution matrix


markerGradient_R_HCC <- data.frame(matrix(nrow = 0, ncol = 3))
markerGradient_NR_HCC <- data.frame(matrix(nrow = 0, ncol = 3))

markerGradient_R_Imm <- data.frame(matrix(nrow = 0, ncol = 3))
markerGradient_NR_Imm <- data.frame(matrix(nrow = 0, ncol = 3))

infiltration_R_Imm <- data.frame(matrix(nrow = 0, ncol = 3))
infiltration_NR_Imm <- data.frame(matrix(nrow = 0, ncol = 3))
infiltration_R_HCC <- data.frame(matrix(nrow = 0, ncol = 3))
infiltration_NR_HCC <- data.frame(matrix(nrow = 0, ncol = 3))


for(col in c('PDL1', 'GranB')){ # remove CCR6, Arg1, CCR6
  
  #col <- 'PDL1'
  NR_data_Imm <- rbind(d2bdry_Imm_7, d2bdry_Imm_8, d2bdry_Imm_9, d2bdry_Imm_6, d2bdry_Imm_27, d2bdry_Imm_28)
  NR_data_Imm$cid <- seq_len(nrow(NR_data_Imm))
  NR_data_HCC <- rbind(d2bdry_HCC_6, d2bdry_HCC_7, d2bdry_HCC_8, d2bdry_HCC_9, d2bdry_HCC_27, d2bdry_HCC_28)
  NR_data_HCC$cid <- seq_len(nrow(NR_data_HCC))
  
  R_data_HCC <- rbind(d2bdry_HCC_22, d2bdry_HCC_23)
  R_data_HCC$cid <- seq_len(nrow(R_data_HCC))
  
  R_data_Imm <- rbind(d2bdry_Imm_22, d2bdry_Imm_23)
  R_data_Imm$cid <- seq_len(nrow(R_data_Imm))
  
  
  
  
  # infiltraing cells
  
  
  
  infiltrating_Imm_R <- R_data_Imm[R_data_Imm$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),] # tumor cells in immune compartment
  infiltrating_Imm_NR <- NR_data_Imm[NR_data_Imm$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),] # tumor cells in immune compartment
  
  infiltrating_HCC_R <- R_data_HCC[R_data_HCC$ctype_no %in% c('ctype1', 'ctype2', 'ctype3', 'ctype5', 'ctype7', 'ctype8', 'ctype10', 'ctype11', 'ctype13', 'ctype15', 'ctype18'),] # immune cells in tumor compartment
  infiltrating_HCC_NR <- NR_data_HCC[NR_data_HCC$ctype_no %in% c('ctype1', 'ctype2', 'ctype3', 'ctype5', 'ctype7', 'ctype8', 'ctype10', 'ctype11', 'ctype13', 'ctype15', 'ctype18'),] # immune cells in tumor compartment
  
  
  
  uninfiltrating_Imm_R <- R_data_Imm[!(R_data_Imm$cid %in% infiltrating_Imm_R$cid),] # tumor cells in immune compartment
  uninfiltrating_Imm_NR <- NR_data_Imm[!(NR_data_Imm$cid %in% infiltrating_Imm_NR$cid),] # tumor cells in immune compartment
  
  uninfiltrating_HCC_R <- R_data_HCC[!(R_data_HCC$cid %in% infiltrating_HCC_R$cid),] # tumor cells in immune compartment
  uninfiltrating_HCC_NR <- NR_data_HCC[!(NR_data_HCC$cid %in% infiltrating_HCC_NR$cid),] # tumor cells in immune compartment
  
  
  
  
  if(col == 'GranB' | col == 'Lag3'){
    R_data_HCC <- R_data_HCC[R_data_HCC$ctype_no == 'ctype8',]
    NR_data_HCC <- NR_data_HCC[NR_data_HCC$ctype_no == 'ctype8',]
    R_data_Imm <- R_data_Imm[R_data_Imm$ctype_no == 'ctype8',]
    NR_data_Imm <- NR_data_Imm[NR_data_Imm$ctype_no == 'ctype8',]
  }
  
  if(col == 'Ki67' | col == 'PDL1' | col == 'CCR6'){
    R_data_HCC <- R_data_HCC[R_data_HCC$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
    NR_data_HCC <- NR_data_HCC[NR_data_HCC$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
    R_data_Imm <- R_data_Imm[R_data_Imm$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
    NR_data_Imm <- NR_data_Imm[NR_data_Imm$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
  }
  
  
  if(col == 'Arg1'){
    R_data_HCC <- R_data_HCC[R_data_HCC$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),]
    NR_data_HCC <- NR_data_HCC[NR_data_HCC$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),]
    R_data_Imm <- R_data_Imm[R_data_Imm$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),]
    NR_data_Imm <- NR_data_Imm[NR_data_Imm$ctype_no %in% c('ctype3', 'ctype5', 'ctype10', 'ctype11', 'ctype2'),]
  }
  
  close_Imm <- data.frame(R_data_Imm[R_data_Imm$dist < 40 , col]) %>%
    mutate(class = 'close-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  far_Imm <- data.frame(R_data_Imm[R_data_Imm$dist >= 40, col]) %>%
    mutate(class = 'far-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  close_HCC <- data.frame(R_data_HCC[R_data_HCC$dist < 40 , col]) %>%
    mutate(class = 'close-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  far_HCC <- data.frame(R_data_HCC[R_data_HCC$dist >= 40, col]) %>%
    mutate(class = 'far-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  
  markerGradient_R_Imm <- rbind(markerGradient_R_Imm, rbind(close_Imm, far_Imm))
  markerGradient_R_HCC <- rbind(markerGradient_R_HCC, rbind(close_HCC, far_HCC))
  
  
  #temp <- data.frame(rbind(close, far))
  
  close_Imm <- data.frame(NR_data_Imm[NR_data_Imm$dist < 40 , col]) %>%
    mutate(class = 'close-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  far_Imm <- data.frame(NR_data_Imm[NR_data_Imm$dist >= 40, col]) %>%
    mutate(class = 'far-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  close_HCC <- data.frame(NR_data_HCC[NR_data_HCC$dist < 40 , col]) %>%
    mutate(class = 'close-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  far_HCC <- data.frame(NR_data_HCC[NR_data_HCC$dist >= 40, col]) %>%
    mutate(class = 'far-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  markerGradient_NR_Imm <- rbind(markerGradient_NR_Imm, rbind(close_Imm, far_Imm))
  markerGradient_NR_HCC <- rbind(markerGradient_NR_HCC, rbind(close_HCC, far_HCC))  
  
  
  #---------------infiltrating cells ----------------------#
  
  
  infil_Imm_R <- data.frame(infiltrating_Imm_R[, col]) %>%
    mutate(class = 'infil-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  uninfil_Imm_R <- data.frame(uninfiltrating_Imm_R[, col]) %>%
    mutate(class = 'uninfil-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  
  infil_Imm_NR <- data.frame(infiltrating_Imm_NR[, col]) %>%
    mutate(class = 'infil-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  uninfil_Imm_NR <- data.frame(uninfiltrating_Imm_NR[, col]) %>%
    mutate(class = 'uninfil-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  
  
  
  infil_HCC_R <- data.frame(infiltrating_HCC_R[, col]) %>%
    mutate(class = 'infil-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  uninfil_HCC_R <- data.frame(uninfiltrating_HCC_R[, col]) %>%
    mutate(class = 'uninfil-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  
  infil_HCC_NR <- data.frame(infiltrating_HCC_NR[, col]) %>%
    mutate(class = 'infil-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  uninfil_HCC_NR <- data.frame(uninfiltrating_HCC_NR[, col]) %>%
    mutate(class = 'uninfil-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
  
  
  infiltration_R_Imm <- rbind(infiltration_R_Imm, rbind(infil_Imm_R, uninfil_Imm_R))
  infiltration_NR_Imm <- rbind(infiltration_NR_Imm, rbind(infil_Imm_NR, uninfil_Imm_NR))
  infiltration_R_HCC <- rbind(infiltration_R_HCC, rbind(infil_HCC_R, uninfil_HCC_R))
  infiltration_NR_HCC <- rbind(infiltration_NR_HCC, rbind(infil_HCC_NR, uninfil_HCC_NR))
  
  
}




markerGradient_R_Imm$response <- 'R'
markerGradient_NR_Imm$response <- 'NR'

markerGradient_R_HCC$response <- 'R'
markerGradient_NR_HCC$response <- 'NR'



markerGradient_Imm <- rbind(markerGradient_R_Imm, markerGradient_NR_Imm)

markerGradient_HCC <- rbind(markerGradient_R_HCC, markerGradient_NR_HCC)

markerGradient <- rbind(markerGradient_Imm, markerGradient_HCC)

stat.test <- markerGradient %>%
  group_by(response, marker) %>%
  wilcox_test(value~class, comparisons = list(c('close-I', 'far-I'), c('close-T', 'far-T'))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()



stat.test <- stat.test %>% add_xy_position(x = "class")

stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '2', 5)
stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '3', 2)
stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '5', 3)

stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '2', 5)
stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '3', 2)
stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '5', 3)

jpeg('D:/DP/Projects/HCC/Figures/Functional_Marker_Profile_Outcome.jpeg', units="in", width=24.5, height=10, res=300)

# adjusted p value

ggboxplot(markerGradient, x = 'class', y = 'value', palette = 'jco', color = 'class', facet.by = c('response', 'marker'), size = 1) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = max(markerGradient$value) + 0.05, label.size = 10) +
  
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 24),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('Normalized Expression') +
  ylim(0, 1.1)

#facet_grid(vars(marker), var(response))
dev.off()

# violin + boxplot



stat.test_prox <- markerGradient %>%
  group_by(class, marker) %>%
  wilcox_test(value~response) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()



stat.test_prox <- stat.test_prox %>% add_xy_position(x = "response")
stat.test_prox

jpeg('D:/DP/Projects/HCC/Figures/Functional_Marker_Profile_Prox.jpeg', units="in", width=10, height=6, res=300)

ggboxplot(markerGradient, x = 'response', y = 'value', palette = 'lancet', color = 'response', 
          bxp.errorbar = TRUE, facet.by = c('marker', 'class')) +
  #stat_boxplot(geom='errorbar', linetype=1, width=0.5)+  #whiskers
  #geom_boxplot(outlier.shape=1)+    
  #theme_bw() +
  #stat_compare_means(aes(label=..p.adj..), comparisons = list(c('close', 'far')), method = 'wilcox.test', label.y = max(markerGradient_R$value) + 0.05,
  #                   p.adjust.methods = 'fdr', size = 5) +
  stat_pvalue_manual(stat.test_prox, label = "p.adj.signif", y.position = max(markerGradient$value) + 0.05, label.size = 6) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 20)) +
  ylab('Normalized Expression') +
  ylim(0, 1.1)
dev.off()


#----- Supplementary Figure: Density plot for LAG3 -------------#

# Imm
DT <- markerGradient_NR_Imm %>%
  filter(class == 'close-I')
DT <- markerGradient_NR_Imm %>%
  filter(class == 'far-I')
DT <- markerGradient_R_Imm %>%
  filter(class == 'close-I')
DT <- markerGradient_R_Imm %>%
  filter(class == 'far-I')

p <- ggplot(DT, aes(x=value, fill= marker)) +
  theme_bw() +
  geom_histogram(position="identity", color='black', alpha = 0.8) +
  scale_fill_manual(values = c('PDL1' = '#cccccc', 'GranB' = '#00cc99')) +
  
  theme(axis.text = element_text(size = 28),
        axis.title = element_blank(),
        legend.position = 'none') +
  ylim(0, 50)
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/NR_far_Imm.png"), width = 8, height = 6, units = "in", dpi = 300)


# HCC
DT <- markerGradient_NR_HCC %>%
  filter(class == 'close-T')
DT <- markerGradient_NR_HCC %>%
  filter(class == 'far-T')
DT <- markerGradient_R_HCC %>%
  filter(class == 'close-T')
DT <- markerGradient_R_HCC %>%
  filter(class == 'far-T')


p <- ggplot(DT, aes(x=value, fill= marker)) +
  theme_bw() +
  geom_histogram(position="identity", color='black', alpha = 0.8) +
  scale_fill_manual(values = c('PDL1' = '#cccccc', 'GranB' = '#00cc99')) +
  theme(axis.text = element_text(size = 28),
        axis.title = element_blank(),
        legend.position = 'none') +
  ylim(0, 250)
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/NR_close_HCC.png"), width = 8, height = 6, units = "in", dpi = 300)



infiltration_R_Imm$response <- 'R'
infiltration_NR_Imm$response <- 'NR'

infiltration_R_HCC$response <- 'R'
infiltration_NR_HCC$response <- 'NR'


infiltration_Imm <- rbind(infiltration_R_Imm, infiltration_NR_Imm)

infiltraiton_HCC <- rbind(infiltration_R_HCC, infiltration_NR_HCC)

infiltration <- rbind(infiltration_Imm, infiltraiton_HCC)

stat.test <- infiltration %>%
  group_by(response, marker) %>%
  wilcox_test(value~class, comparisons = list(c('infil-I', 'uninfil-I'), c('infil-T', 'uninfil-T'))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()



stat.test <- stat.test %>% add_xy_position(x = "class")

stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '2', 5)
stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '3', 2)
stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '5', 3)

stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '2', 5)
stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '3', 2)
stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '5', 3)

jpeg('D:/DP/Projects/HCC/Figures/Functional_Marker_Infiltration.jpeg', units="in", width=24.5, height=10, res=300)

# adjusted p value

ggboxplot(infiltration, x = 'class', y = 'value', palette = 'jco', color = 'class', facet.by = c('response', 'marker'), size = 1) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = max(markerGradient$value) + 0.05, label.size = 10) +
  
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 24),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('Normalized Expression') +
  ylim(0, 1.1)

#facet_grid(vars(marker), var(response))
dev.off()














# ------------------------------------------------------------ #
#           Boxplot - Cell Type - Close versus Far             #
# ------------------------------------------------------------ #



stat.test <- simple_Imm_in_Tumor_all %>%
  wilcox_test(density~response) %>%
  add_significance()

p <- ggplot(simple_Imm_in_Tumor_all, aes(x = response, y = density), fill = 'black', fill = 'transparent') +
  geom_violin(position = 'dodge') +
  theme_bw() +
  geom_jitter(shape = 21, size = 6, aes(fill = response), alpha = 0.4)+
  #stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(Imm_in_Tumor_all$Freq) + 8, label.size = 10) +
  geom_boxplot(aes(ymin = min(density),ymax = max(density)), fill = 'transparent', outlier.size = -1) +
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 26, angle = 90),
        axis.title.x = element_blank(),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white', size = 2),
        #panel.spacing = unit(0, "mm"),
        axis.line = element_blank()) +
  scale_fill_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  ylab("") 
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Figure4G_1.png"), width = 6, height = 10, units = "in", dpi = 300)



stat.test <- simple_Tumor_in_Imm_all %>%
  wilcox_test(density~response) %>%
  add_significance()

p <- ggplot(simple_Tumor_in_Imm_all, aes(x = response, y = density), fill = 'black', fill = 'transparent') +
  geom_violin(position = 'dodge') +
  theme_bw() +
  geom_jitter(shape = 21, size = 6, aes(fill = response), alpha = 0.4)+
  #stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(Tumor_in_Imm_all$Freq) + 8, label.size = 10) +
  geom_boxplot(aes(ymin = min(density),ymax = max(density)), fill = 'transparent', outlier.size = -1) +
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 26, angle = 90),
        axis.title.x = element_blank(),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white', size = 2),
        #panel.spacing = unit(0, "mm"),
        axis.line = element_blank()) +
  scale_fill_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  ylab("")
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Figure4G_2.png"), width = 6, height = 10, units = "in", dpi = 300)



#----------- Just bar plot ---------------#

p <- ggbarplot(simple_Imm_in_Tumor_all, x = "core", y = "density",
               color = "core", fill = 'transparent', orientation = 'horiz',
               palette = c('#21b7bd', '#21b7bd', '#ec776d', '#ec776d', '#ec776d', '#ec776d', '#ec776d', '#ec776d'),
               order = c('22', '23', '6', '7', '8', '9', '27', '28'),
               size = 1.5, position = position_dodge()) +
  rremove("legend") +
  font("xy.text", size = 24) +
  font('xy.title', size = 26) +
  xlab('Core') +
  ylab("")

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Figure4G_1.png"), width = 10, height = 6, units = "in", dpi = 300)






p <- ggbarplot(simple_Tumor_in_Imm_all, x = "core", y = "density",
          color = "core", fill = 'transparent', orientation = 'horiz',
          palette = c('#21b7bd', '#21b7bd', '#ec776d', '#ec776d', '#ec776d', '#ec776d', '#ec776d', '#ec776d'),
          order = c('22', '23', '6', '7', '8', '9', '27', '28'),
          size = 1.5, position = position_dodge(), xlab = FALSE) +
  rremove("legend") +
  rremove("x.title") +
  font("xy.text", size = 24) +
  font('xy.title', size = 26) +
  xlab('Response') +
  ylab("")
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Figure4G_2.png"), width = 10, height = 6, units = "in", dpi = 300)























# ----------------------------------------------- #
#           Bar plot  - Supplementary             #
# ----------------------------------------------- #

{

Core_dat <- data.frame(matrix(nrow = 4, ncol = 3))

core <- 6
d2bdry_Imm <- eval(sym(paste0('d2bdry_Imm_',core)))
count_Imm_close <- nrow(d2bdry_Imm[d2bdry_Imm$dist < 40,])
count_Imm_far <- nrow(d2bdry_Imm[d2bdry_Imm$dist >= 40,])

d2bdry_HCC <- eval(sym(paste0('d2bdry_HCC_',core)))
count_HCC_close <- nrow(d2bdry_HCC[d2bdry_HCC$dist < 40,])
count_HCC_far <- nrow(d2bdry_HCC[d2bdry_HCC$dist >= 40,])


# ------- combine DF -------#

Core_dat[1,] <- c('Immune', 'Close', count_Imm_close)
Core_dat[2,] <- c('Immune', 'Far', count_Imm_far)
Core_dat[3,] <- c('Tumor', 'Close', count_HCC_close)
Core_dat[4,] <- c('Tumor', 'Far', count_HCC_far)
colnames(Core_dat) <- c('compartment', 'RR', 'value')
Core_dat$value <- as.numeric(Core_dat$value)

bp <- ggplot(Core_dat, aes(x = RR, y = value, color = compartment)) +   
  theme_bw() +
  geom_bar(stat = 'identity', width = 0.7, position = position_dodge(width = 0.8), fill = 'white', size = 2) +
  scale_color_manual(values = c('Tumor' = 'red', 'Immune' = 'blue')) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 32),
        axis.title = element_text(size = 34),
        legend.position = 'none') +
  ylim(0,500) +
  ylab('Cell counts')
bp
ggsave(bp, file= paste0('D:/DP/Projects/HCC/Figures/compositon_bp_', core, '.jpeg'), width = 8, height = 7, units = "in", dpi = 300)
}
#------------- Box plot ---------------#
{
  markerGradient_HCC <- data.frame(matrix(nrow = 0, ncol = 3))
  markerGradient_Imm <- data.frame(matrix(nrow = 0, ncol = 3))
  
  for(col in c('CCR6', 'Arg1', 'PDL1', 'HLADR', 'GranB', 'Ki67')){
    
    #col <- 'PDL1'
    data_Imm <- d2bdry_Imm_6
    
    data_HCC <- d2bdry_HCC_6
    
    
    if(col == 'GranB'){
      data_HCC <- data_HCC[data_HCC$ctype_no == 'ctype8',]
      
      data_Imm <- data_Imm[data_Imm$ctype_no == 'ctype8',]
    }
    
    if(col == 'Ki67' | col == 'PDL1' | col == 'CCR6'){
      data_HCC <- data_HCC[data_HCC$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
      data_Imm <- data_Imm[data_Imm$ctype_no %in% c('ctype9', 'ctype14', 'ctype17', 'ctype12', 'ctype16'),]
    }
    
    
    
    close_Imm <- data.frame(data_Imm[data_Imm$dist < 40 , col]) %>%
      mutate(class = 'close-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
    far_Imm <- data.frame(data_Imm[data_Imm$dist >= 40, col]) %>%
      mutate(class = 'far-I', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
    
    close_HCC <- data.frame(data_HCC[data_HCC$dist < 40 , col]) %>%
      mutate(class = 'close-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
    far_HCC <- data.frame(data_HCC[data_HCC$dist >= 40, col]) %>%
      mutate(class = 'far-T', marker = col) %>% `colnames<-`(c("value", "class", "marker"))
    
    
    markerGradient_Imm <- rbind(markerGradient_Imm, rbind(close_Imm, far_Imm))
    markerGradient_HCC <- rbind(markerGradient_HCC, rbind(close_HCC, far_HCC))
    
    
  }
  
  
  markerGradient <- rbind(markerGradient_Imm, markerGradient_HCC)
  
  stat.test <- markerGradient %>%
    group_by(marker) %>%
    wilcox_test(value~class, comparisons = list(c('close-I', 'far-I'), c('close-T', 'far-T'))) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance()
  
  
  
  stat.test <- stat.test %>% add_xy_position(x = "class")
  
  stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '2', 5)
  stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '3', 2)
  stat.test$xmin <- replace(stat.test$xmin, stat.test$xmin == '5', 3)
  
  stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '2', 5)
  stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '3', 2)
  stat.test$xmax <- replace(stat.test$xmax, stat.test$xmax == '5', 3)
  
  
  
  
  jpeg('D:/DP/Projects/HCC/Figures/Core6_Boxplot.jpeg', units="in", width=18, height=6, res=300)
  
  ggboxplot(markerGradient, x = 'class', y = 'value', palette = 'jco', color = 'class', facet.by = c('marker')) +
    #stat_boxplot(geom='errorbar', linetype=1, width=0.5)+  #whiskers
    #geom_boxplot(outlier.shape=1)+    
    #theme_bw() +
    #stat_compare_means(aes(label=..p.adj..), comparisons = list(c('close', 'far')), method = 'wilcox.test', label.y = max(markerGradient_R$value) + 0.05,
    #                   p.adjust.methods = 'fdr', size = 5) +
    stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = max(markerGradient$value) + 0.05, label.size = 10) +
    theme(axis.title = element_text(size = 26),
          axis.text = element_text(size = 24),
          axis.title.x = element_blank(),
          legend.position = 'none',
          strip.text = element_text(size = 24),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(cols = vars(marker)) +
    ylab('Normalized Expression') +
    ylim(0, 1.1)
  dev.off()
  
  
  
}
