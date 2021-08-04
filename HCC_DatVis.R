# Hepatocelluler Carcinoma SpatStat Script #
# @Author: Haoyang Mi


library(ggplot2)
library(spatstat)
library(umap)
library(rmarkdown)
library(Rtsne)
library(dplyr)
library(FlowSOM)
library(flowCore)
library(flowAssist)
library(plotly)
library(ggvoronoi)
library(matrixStats)

setwd("D:/DP/Projects/HCC")

# read data

HCCdata <- readRDS('D:/DP/Data/HCC/hccdataset')
unique(HCCdata$ctype)
#remove UA Noncell




test <- HCCdata %>%
  filter(ctype != 'UA Noncell') %>%
  filter(response == 'NR') %>%
  group_by(Core, ctype) %>%
  tally() %>%
  mutate(density = n / (pi*(0.4)^2))

for(id in unique(test$ctype)) {
  #id <- 'Imm DP T'
  dat <- test[test$ctype == id,]
  
  print(paste('#----------', id))
  print(paste('Max is:', max(dat$density)))
  print(paste('min is:', min(dat$density)))
  print(paste('Mean is:', mean(dat$density)))
  print(paste('SD is:', sd(dat$density)))
}
#-------------------------------------#
#                 tSNE                #
#-------------------------------------#




# gradient color - marker label

#cofactor <- 5
#HCCdata <- asinh(HCCdata[,1:32]/cofactor)


# scale using quantile
rng <- colQuantiles(as.matrix(HCCdata[,1:32]), probs = c(0.01, 0.99))


HCCdata0 <- t((t(HCCdata[,1:32]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
HCCdata0[HCCdata0 < 0] <- 0; HCCdata0[HCCdata0 > 1] <- 1
HCCdata0  <- data.frame(HCCdata0)


for(marker in colnames(HCCdata_int[,1:27])){
  
  tempplot <- ggplot(data = tSNE_all.layout, aes(X1, X2, color = HCCdata0[,marker]))+
    geom_point(size = 1) +
    theme_classic() +
    scale_color_gradientn(colours = c('#263e9d', '#dbecfa', '#a00d2f')) +
    theme(legend.position = 'none',
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_line(size = 2),
          aspect.ratio = 1) 
  
  tempplot
  ggsave(tempplot, file = paste0('./Figures/tSNE-allMarkers/tSNE_', marker,'.jpeg'), width = 10, height = 10, units = "in", dpi = 300)
  
}


# UMAP for Responder
HCCdata0$response <- HCCdata$response
HCCdata0$ctype <- HCCdata$ctype
HCCdata0$ctype_no <- HCCdata$ctype_no

HCCdata0 <- HCCdata0[HCCdata0$ctype != 'UA Noncell',]

HCCdata_int <- select(HCCdata0,'HisH3', 'DNA', 'HLADR', 'aSMA', 'Collagen',
                      'Vimentin', 'Ecad', 'PanK', 'CCR6', 'Arg1', 'CD45RO',
                      'CD3', 'CD45', 'CD4', 'GranB', 'CD28', 'PDL1', 'Foxp3', 'Lag3', 
                      'Casp3', 'CD8a', 'CD20', 'CD68', 'CD16', 'CD163', 'CD15', 'Ki67', 'response')

HCCdata_int_R <- HCCdata_int[HCCdata_int$response == 'R',1:27]
HCCdata_int_NR <- HCCdata_int[HCCdata_int$response == 'NR',1:27]

tSNE_all <- Rtsne(HCCdata_int, check_duplicates = FALSE, pca = TRUE, perplexity = 100, dims = 2, max_iter = 5000)
tSNE_all.layout <- data.frame(tSNE_all$Y)
all.label <- HCCdata0$ctype_no

  
jpeg('./Figures/tSNE_all.jpeg', unit = 'in', width = 10, height = 10, res = 300)
ggplot(data = tSNE_all.layout, aes(X1, X2, color = as.factor(all.label)))+
  geom_point(size = 1) +
  theme_void() +
  scale_color_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        aspect.ratio = 1) 
dev.off()

tSNE_R <- Rtsne(HCCdata_int_R, check_duplicates = FALSE, pca = TRUE, perplexity = 100, dims = 2, max_iter = 5000)
tSNE_R.layout <- data.frame(tSNE_R$Y)
R.label <- HCCdata0[HCCdata0$response == 'R',]$ctype_no


tSNE_NR <- Rtsne(HCCdata_int_NR, check_duplicates = FALSE, pca = TRUE, perplexity = 100, dims = 2, max_iter = 5000)
tSNE_NR.layout <- data.frame(tSNE_NR$Y)

NR.label <- HCCdata0[HCCdata0$response == 'NR',]$ctype_no

jpeg('./Figures/tSNE_R.jpeg', unit = 'in', width = 10, height = 10, res = 300)
ggplot(data = tSNE_R.layout, aes(X2, X1, color = as.factor(R.label)))+
  geom_point(size = 1) +
  theme_void() +
  scale_color_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        aspect.ratio = 1) 
  

dev.off()

jpeg('./Figures/tSNE_NR.jpeg', unit = 'in', width = 10, height = 10, res = 300)
ggplot(data = tSNE_NR.layout, aes(X1, X2, color = as.factor(NR.label)))+
  geom_point(size = 1) +
  theme_void() +
  scale_color_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        aspect.ratio = 1) 
dev.off()


#-----------------------------------#
#               UMAP                #
#-----------------------------------#




R.umap = umap(HCCdata_int_R)

plot(R.umap$layout)


#---------------------------------#
#             FlowSOM             #
#---------------------------------#


FlowSOM_dat <- DFtoFF(HCCdata_int[, 1:31])
fSOM <- FlowSOM(FlowSOM_dat,
                  # Input options:
                  #compensate = TRUE, transform = TRUE, toTransform=c(1:31),
                  scale = TRUE,
                  # SOM options:
                  colsToUse = c(1:31), xdim = 7, ydim = 7,
                  # Metaclustering options:
                   nClus = 18)
PlotStars(fSOM[[1]],
             backgroundValues = as.factor(fSOM[[2]]))


#-----------------------------------------#
#       cell type-level characterization       #
#-----------------------------------------#
pChange_all <- data.frame(matrix(nrow = 0, ncol = 0))
for (ctype in unique(HCCdata$ctype_no)) {
  
  
  #ctype <- 'ctype15'
  
  rHCC <- HCCdata[HCCdata$ctype_no == ctype & HCCdata$response == 'R',]
  
  nrHCC <- HCCdata[HCCdata$ctype_no == ctype & HCCdata$response == 'NR',]
  
  ctype_name <- as.character(unique(rHCC$ctype))
  cellCount_vec <- vector()
  
  for (core in unique(rHCC$Core)) {
    
    cellCount <- nrow(rHCC[rHCC$Core == core,])
    
    cellCount_vec <- c(cellCount_vec, cellCount)
  }
  
  meanR <- mean(cellCount_vec)
  sdR <- sd(cellCount_vec)
  
  cellCount_vec <- vector()
  
  for (core in unique(nrHCC$Core)) {
    
    cellCount <- nrow(nrHCC[nrHCC$Core == core,])
    
    cellCount_vec <- c(cellCount_vec, cellCount)
  }
  meanNR <- mean(cellCount_vec)
  sdNR <- sd(cellCount_vec)
  
  pChange <- (meanR - meanNR)*100/meanNR
  
  pChange_all <- rbind(pChange_all, cbind(ctype_name, as.numeric(pChange)))
  
}

write.csv(pChange_all, './Data/pChange_all.csv', row.names = FALSE)


#-----------------------------------------#
#       core-level characterization       #
#-----------------------------------------#


ncore <- max(HCCdata$Core)

nctype <- length(unique(HCCdata$ctype))

allCoreDat <- data.frame(matrix(nrow = 0, ncol = 0))

for (core_id in 1:ncore) {

  dens_vec <- data.frame(matrix(nrow = 1, ncol = 0))
  per_vec <- data.frame(matrix(nrow = 1, ncol = 0))
  
  #core_id <- 1
  subDat <- HCCdata[HCCdata$Core == core_id,] 
    
  
  for(type_id in 1:nctype){
    
    #type_id <- 1
    
    type_id <- paste('ctype', type_id, sep = '')
    
    ssDat <- subDat[subDat$ctype_no == type_id,]
    
    # calculate cell type density
    dens <- nrow(ssDat)/(pi*0.75^2)
    percent <- nrow(ssDat)/nrow(subDat)
    
    # density vector
    dens_vec <- cbind(dens_vec, dens)
    per_vec <- cbind(per_vec, percent)
    

  }
 
  core_vec <- cbind(core_id, dens_vec, unique(subDat$response)) 
  
  percent_core_vec <- cbind(core_id, dens_vec, unique(subDat$response)) 

  allCoreDat <- rbind(allCoreDat, core_vec)
  allCoreDat_percent <- rbind(allCoreDat, percent_core_vec)
  
}


colnames(allCoreDat) <- c('Core', paste('ctype',seq(1,18),'_dens', sep = ''), 'response')

colnames(allCoreDat_percent) <- c('Core', paste('ctype',seq(1,18),'_percent', sep = ''), 'response')


write.csv(allCoreDat, 'D:/DP/Data/HCC/allCoreDat.csv')
write.csv(allCoreDat_percent, 'D:/DP/Data/HCC/allCoreDat_percent.csv')

#--------------------------------------#
#               Finish                 #
#--------------------------------------#


#-------------------------------------------#
#           assign patient to core          #
#-------------------------------------------#

allCoreDat$Patient <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3), rep(9,4), rep(10,3), rep(11,3), rep(12,3))


# patient-level merge

ptDF <- data.frame(matrix(nrow = 0, ncol = 0))
for(p in 1:12){
  
  
  # retrive patient dat
  
  ptData <- allCoreDat[allCoreDat$Patient == p,]
  
  response <- as.character(unique(ptData$response))
  # numerical estimates average
  ptVec <- t(colMeans(ptData[,2:19]))
  
  
  ptDF <- rbind(ptDF, cbind(p, ptVec, response))
}




#------------------------------#
#         Core Data Heatmap    #
#------------------------------#


Core_Heatmap <- data.frame(matrix(nrow = 0, ncol = 18))
Vec <- rep(0,17)
Name <- c('ctype1', 'ctype2', 'ctype3', 'ctype5', 'ctype6','ctype7',
           'ctype8', 'ctype9', 'ctype10', 'ctype11', 'ctype12','ctype13',
           'ctype14', 'ctype15', 'ctype16', 'ctype17', 'ctype18')
names(Vec) <- Name


for(core in seq_len(37)){
  
  #core <- 4
  Dat <- HCCdata[HCCdata$Core == core,]
  
  Dat_Table <- t(data.frame(table(Dat$ctype_no)))
  
  
  Dat_Table_name <- as.character(Dat_Table[1,])
  
  insert <- t(as.numeric(as.character(Dat_Table[2,])))
  colnames(insert) <- Dat_Table_name
  
  # join two vectors by column names
  datC <- full_join(data.frame(t(Vec)), data.frame(insert), copy = FALSE)
  
  Core_Heatmap <- rbind(Core_Heatmap, cbind(core, datC[2,]))
  
}

Core_Heatmap[is.na(Core_Heatmap)] <- 0

write.csv(Core_Heatmap, './Data/Core_Heatmap.csv', row.names = FALSE)




#------------------------------#
#           Point Patterns     #
#------------------------------#

# Get R and NR dat

Sample_dat <- HCCdata[HCCdata$Core == 22,]


jpeg('./Figures/Core_22_PP_phenotyped.jpeg', unit = 'in', width = 5, height = 5, res = 300)
ggplot(data = Sample_dat, aes(Xcoord/1000, Ycoord/1000, fill = as.factor(ctype_no)))+
  geom_point(size = 6, shape = 21) +
  theme_bw() +
  theme(legend.position = 'none',
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    #axis.line = element_line(size = 1),
    aspect.ratio = 1) +
  scale_fill_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                                'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                                'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                                'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                                'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                                'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  xlab('X, mm') +
  ylab('Y, mm')
dev.off()


jpeg('./Figures/Core_22_PP_phenotyped_sub.jpeg', unit = 'in', width = 5, height = 5, res = 300)
ggplot(data = Sample_dat, aes(Xcoord/1000, Ycoord/1000, fill = as.factor(ctype_no)))+
  geom_point(size = 6, shape = 21) +
  theme_bw() +
  xlim(0.38,0.58) +
  ylim(0.26, 0.46) +
  theme(legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
        aspect.ratio = 1) +
  scale_fill_manual(values = c('ctype1' = '#173ddd', 'ctype2' = '#fcc5f5', 'ctype3' = '#b346ce',
                               'ctype5' = '#86208e', 'ctype6' = '#83e276', 'ctype7' = '#6df2f2',
                               'ctype8' = '#94cfea', 'ctype9' = '#ed8a48', 'ctype10' = '#d382f9',
                               'ctype11' = '#ddaaf7', 'ctype12' = '#000000', 'ctype13' = '#e019ba',
                               'ctype14' = '#f43a1c', 'ctype15' = '#3573d8', 'ctype16' = '#b7b7b7',
                               'ctype17' = '#f97862', 'ctype18' = '#5d95e8')) +
  xlab('X, mm') +
  ylab('Y, mm')
dev.off()
#-------------------------------------#
#           Basic comparison          #
#-------------------------------------#

# Get R and NR dat

HCC_R <- ptDF[ptDF$response == 'R',]

HCC_NR <- ptDF[ptDF$response == 'NR',]

write.csv(HCC_R, 'allCoreDat_R.csv')
write.csv(HCC_NR, 'allCoreDat_NR.csv')

# p.value and multiple comparison correction

pval <- data.frame(matrix(nrow = 0, ncol = 0))

for(col in 2:19){
  res <-  wilcox.test(as.numeric(as.character(HCC_R[,col])), as.numeric(as.character(HCC_NR[,col])))
  pval <- rbind(pval, res$p.value)
  print(res$p.value)
}

p.adjust(as.numeric(pval[,1]), method = 'BY')



#----------------------------------------------#
#           Cor - Imm vs Hep vs Macro          #
#----------------------------------------------#

count <- data.frame(matrix(nrow = 0, ncol = 6))
for(core in seq_len(37)){
  
  #core <- 1
  Dat <- HCCdata[HCCdata$Core == core,]
  
  cImm <- nrow(Dat[paste('ctype', c(1, 2, 3, 5, 7, 8, 10, 11, 13, 15, 18), sep = '') %in% Dat$ctype_no,])
  
  cHCC <- nrow(Dat[paste('ctype', c(9, 12, 14, 16, 17), sep = '') %in% Dat$ctype_no,])
  
  cLym <- nrow(Dat[paste('ctype', c(1, 7, 8, 15, 18), sep = '') %in% Dat$ctype_no, ])
  
  cMye <- nrow(Dat[paste('ctype', c(2, 3, 5, 10, 11, 13), sep ='') %in% Dat$ctype_no, ])
  
  Response <- as.character(unique(Dat$response))
  
  count <- rbind(count, cbind(core, cImm, cHCC, cLym, cMye, Response))
}


R <- count[count$Response == 'R',]

NR <- count[count$Response == 'NR',]

count$cImm <- as.numeric(as.character(count$cImm))
count$cHCC <- as.numeric(as.character(count$cHCC))
count$cLym <- as.numeric(as.character(count$cLym))
count$cMye <- as.numeric(as.character(count$cMye))

count_R <- count[count$Response == 'R',]
count_NR <- count[count$Response == 'NR',]

cor(count_R$cMye, count_R$cHCC, method = 'spearman')
cor(count_NR$cMye, count_NR$cHCC, method = 'spearman')

jpeg('./Figures/Cor-Imm-Hep.jpeg', unit = 'in', width = 5, height = 5, res = 300)
ggplot(data = count, aes(cHCC, cImm, fill = as.factor(Response)))+
  stat_smooth(aes(color = as.factor(Response)), method='lm', se = FALSE, size = 2) +
  geom_point(size = 8, shape = 21) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1),
    aspect.ratio = 1) +
  scale_color_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  scale_fill_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  xlab('HCC counts') +
  ylab('Immune cell counts')
  
dev.off()


jpeg('./Figures/Cor-Lym-Mye.jpeg', unit = 'in', width = 5, height = 5, res = 300)
ggplot(data = count, aes(cLym, cMye, fill = as.factor(Response)))+
  stat_smooth(aes(color = as.factor(Response)), method='lm', se = FALSE, size = 2) +
  geom_point(size = 8, shape = 21) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1),
    aspect.ratio = 1) +
  scale_color_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  scale_fill_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  xlab('Lymphoid cell counts') +
  ylab('Myeloid cell counts')

dev.off()

jpeg('./Figures/Cor-HCC-Lym.jpeg', unit = 'in', width = 5, height = 5, res = 300)
ggplot(data = count, aes(cHCC, cLym, fill = as.factor(Response)))+
  theme_bw() +
  stat_smooth(aes(color = as.factor(Response)), method='lm', se = FALSE, size = 2) +
  geom_point(size = 8, shape = 21) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1),
    aspect.ratio = 1) +
  scale_color_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  scale_fill_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  xlab('HCC counts') +
  ylab('Lymphoid cell counts')

dev.off()

cor(R$cHCC, R$cMye, method = 'spearman')
cor(NR$cHCC, NR$cMye, method = 'spearman')

jpeg('./Figures/Cor-HCC-Mye.jpeg', unit = 'in', width = 5, height = 5, res = 300)
ggplot(data = count, aes(cHCC, cMye, fill = as.factor(Response)))+
  theme_bw() +
  stat_smooth(aes(color = as.factor(Response)), method='lm', se = FALSE, size = 2) +
  geom_point(size = 8, shape = 21) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 1),
    aspect.ratio = 1) +
  scale_color_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  scale_fill_manual(values = c('R' = '#21b7bd', 'NR' = '#ef776d')) +
  xlab('HCC counts') +
  ylab('Myeloid cell counts')

dev.off()



#---------------------------#
#          Pie Chart        #
#---------------------------#
PieData <- as.data.frame(table(HCCdata$ctype_no))
PieData$Perc <- round(PieData$Freq*100/sum(PieData$Freq), digit = 1)


p <- plot_ly(PieData, labels = ~Var1, values = ~Perc, text = ~paste(Perc, '%', sep = ''), type = 'pie',textposition = 'outside',
             textinfo = 'text',
             textfont = list(size = 40, family = 'Arial'),
             marker = list(colors = c('#173ddd', '#d382f9' , '#ddaaf7', '#000000', '#e019ba',
                                      '#f43a1c', '#3573d8', '#b7b7b7', '#f97862', '#5d95e8',
                                      '#fcc5f5', 'ctype3' = '#b346ce', '#86208e', '#83e276','#6df2f2', '#94cfea', '#ed8a48' ))) %>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend = FALSE,
         autosize = F)

p

orca(p, './Figures/PoW-All.svg')

# responders
HCCdataR <- HCCdata[HCCdata$response == 'R', ]

PieDataR <- as.data.frame(table(HCCdataR$ctype_no))
PieDataR$Perc <- round(PieDataR$Freq*100/sum(PieDataR$Freq), digit = 1)


p <- plot_ly(PieDataR, labels = ~Var1, values = ~Perc, text = ~paste(Perc, '%', sep = ''), type = 'pie',textposition = 'outside',
             textinfo = 'text',
             textfont = list(size = 40, family = 'Arial'),
             marker = list(colors = c('#173ddd', '#d382f9' , '#ddaaf7', '#000000', '#e019ba',
                                      '#f43a1c', '#3573d8', '#b7b7b7', '#f97862', '#5d95e8',
                                      '#fcc5f5', 'ctype3' = '#b346ce', '#86208e', '#83e276','#6df2f2', '#94cfea', '#ed8a48' ))) %>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend = FALSE,
         autosize = F)

p

orca(p, './Figures/PoW-R.svg')

# Non-Responders

HCCdataNR <- HCCdata[HCCdata$response == 'NR', ]

PieDataNR <- as.data.frame(table(HCCdataNR$ctype_no))

PieDataNR$Perc <- round(PieDataNR$Freq*100/sum(PieDataNR$Freq), digit = 1)


p <- plot_ly(PieDataNR, labels = ~Var1, values = ~Perc, text = ~paste(Perc, '%', sep = ''), type = 'pie',textposition = 'outside',
             textinfo = 'text',
             textfont = list(size = 40, family = 'Arial'),
             marker = list(colors = c('#173ddd', '#d382f9' , '#ddaaf7', '#000000', '#e019ba',
                                      '#f43a1c', '#3573d8', '#b7b7b7', '#f97862', '#5d95e8',
                                      '#fcc5f5', 'ctype3' = '#b346ce', '#86208e', '#83e276','#6df2f2', '#94cfea', '#ed8a48' ))) %>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend = FALSE,
         autosize = F)

p

orca(p, './Figures/PoW-NR.svg')
#----------------------------#
#          Histograms        #
#----------------------------#

jpeg('./Figures/HCC-Barplot.jpeg', unit = 'in', width = 8, height = 5, res = 300)
ggplot(data = count, aes(x = core, y = cHCC, width = 0.8))+
  theme_classic() +
  geom_bar(stat = 'identity', color = 'black', fill = '#cfdae2') +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    aspect.ratio = 1)
  
dev.off()

jpeg('./Figures/Mye-Barplot.jpeg', unit = 'in', width = 8, height = 5, res = 300)
ggplot(data = count, aes(x = core, y = cMye, width = 0.8))+
  theme_classic() +
  geom_bar(stat = 'identity', color = 'black', fill = '#cfdae2') +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    aspect.ratio = 1)

dev.off()

jpeg('./Figures/Lym-Barplot.jpeg', unit = 'in', width = 8, height = 5, res = 300)
ggplot(data = count, aes(x = core, y = cLym, width = 0.8))+
  theme_classic() +
  geom_bar(stat = 'identity', color = 'black', fill = '#cfdae2') +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    aspect.ratio = 1)

dev.off()


jpeg('./Figures/Imm-Barplot.jpeg', unit = 'in', width = 8, height = 5, res = 300)
ggplot(data = count, aes(x = core, y = cImm, width = 0.8))+
  theme_classic() +
  geom_bar(stat = 'identity', color = 'black', fill = '#cfdae2') +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    aspect.ratio = 1)

dev.off()


#--------------------------------------------------------------#
#           Pairwise comparison in R and NR - Heatmap          #
#--------------------------------------------------------------#


# pvalue vector
pvalue <- c()
# output from CoreHeatmap
for(col in seq_len(17)){
  for(row in seq_len(17)){
    
    #col <- 1
    row <- 2
    # col and row cell type
    colType <- paste('ctype', col, sep = '')
    rowType <- paste('ctype', row, sep = '')
    
    
    colVec <- Core_Heatmap[colType]
    rowVec <- Core_Heatmap[rowType]
    
    wilcox.test(colVec[,1], rowVec[,1])

    
    
  }
}







#------------------------------------#
#           Graph Embedding          #
#------------------------------------#
GNN_vec <- read.csv('./Data/Features/2/nci.csv', row.names = 1)

tSNE_GNN <- Rtsne(GNN_vec, check_duplicates = FALSE, pca = TRUE, perplexity = 1, dims = 2, max_iter = 1000)

tSNE_GNN.layout <- data.frame(tSNE_GNN$Y)
all.label <- seq_len(nrow(GNN_vec))


jpeg('./Figures/tSNE_all.jpeg', unit = 'in', width = 10, height = 10, res = 300)
ggplot(data = tSNE_GNN.layout, aes(X1, X2, color = as.factor(all.label)))+
  geom_point(size = 8) +
  theme_classic() +
  theme(#legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(size = 1),
        aspect.ratio = 1) +
  geom_text(data = tSNE_GNN.layout, aes(x = X1, y = X2, label = seq_len(26)), color = 'black', size = 6) 
dev.off()



