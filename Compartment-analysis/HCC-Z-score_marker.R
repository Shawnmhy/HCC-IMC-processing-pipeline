# Hepatocelluler Carcinoma SpatStat Script #
# @Author: Haoyang Mi
# Z score
library(rmarkdown)
library(spatstat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RANN)
library(matrixStats)
library(radiomics)
library(factoextra)
library(FactoMineR)

setwd("D:/DP/Projects/HCC")


# read functions
source('D:/DP/Projects/HCC/Functions.r')
# read data

HCCdata <- readRDS("D:/DP/Data/HCC/hccdataset")

unique(HCCdata$ctype)
#remove UA Noncell

# Include certain markers
HCCdata_ <- HCCdata[,c('CD45', 'CD16', 'CD163', 'CD68', 'CD3', 'CD4', 'CD8a', 'CD20', 'CD28', 'Foxp3', 'CD45RO', 'CD15',
                       'CCR6', 'Arg1', 'PanK', 'Ecad', 'PDL1', 'Ki67', 'HLADR', 'GranB', 'Casp3', 'Lag3', 'Collagen', 'Vimentin',
                       'Core', 'Xcoord', 'Ycoord')]


#temp <- HCCdata[order(HCCdata$HisH3, decreasing = TRUE),]
#temp[round(0.9*n), 'HisH3']

mType <- c('CD45', 'CD16', 'CD163', 'CD68', 'CD3', 'CD4', 'CD8a', 'CD20', 'CD28', 'Foxp3', 'CD45RO', 'CD15',
           'CCR6', 'Arg1', 'PanK', 'Ecad', 'PDL1', 'Ki67', 'HLADR', 'GranB', 'Casp3', 'Lag3', 'Collagen', 'Vimentin')


#FlowSOM_dat <- DFtoFF(HCCdata_[, 1:24])

#temp <- fSOM$data
#mean(HCCdata[HCCdata$ctype_no == 'ctype4', 'DNA'])
# arcsinh transformation

#cofactor <- 5
#HCCdata_[,1:24] <- asinh(HCCdata_[,1:24]/cofactor)
HCCdata_$ctype <- HCCdata$ctype
HCCdata_$ctype_no <- HCCdata$ctype_no

# scale using quantile
rng <- colQuantiles(as.matrix(HCCdata_[,1:24]), probs = c(0.01, 0.99))


HCCdata0 <- t((t(HCCdata_[,1:24]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
HCCdata0[HCCdata0 < 0] <- 0; HCCdata0[HCCdata0 > 1] <- 1

HCCdata_[,1:24] <- HCCdata0
HCCdata_ <- HCCdata_ %>%
  filter(ctype != 'UA Noncell')
#HCCdata_[,1:24] <- apply(HCCdata_[,1:24], 2, function(x){ (x - min(x)) / (max(x) - min(x)) })


#hist(as.numeric(HCCdata_$Lag3))


textureFeature_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(core in seq_len(37)){
  
  # search neighbors
  core <- 8
  coreDat <- HCCdata_[HCCdata_$Core == core,]
  
  pos <- coreDat[,c('Xcoord', 'Ycoord')]
  dist <- nn2(pos, query = pos, k = nrow(pos), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  nn.idx <- data.frame(dist$nn.idx)
  
  nn.idx$source <- seq_len(nrow(coreDat))
  
  nn.merge <- melt(nn.idx, id.vars = 'source')
  
  # clear non-valid rows
  nn.merge <- nn.merge[nn.merge$value != 0, -2]
  
  
  # remove itself
  
  #nn.merge <- nn.merge %>%
  #  mutate(source = pmin(source, value),
  #         target = pmax(source, value)) %>%
  #  distinct(source, target)
  
  nn.merge$diff <- nn.merge$source - nn.merge$value
  nn.merge <- nn.merge[nn.merge$diff != 0, -3]
  
  
  # Create and fill heatmap
  Heatmap <- data.frame(matrix(nrow = 24, ncol = 24))
  
  colnames(Heatmap) <- mType
  rownames(Heatmap) <- mType
  
  
  grid <- expand.grid(mType, mType)
  for(id in seq_len(nrow(grid))){
    print(id)
    sType <- as.character(grid[id, 1])
    tType <- as.character(grid[id, 2])
    
    
    nn.merge.val <- nn.merge
    
    nn.merge.val$source <- coreDat[nn.merge$source,sType]
    nn.merge.val$value <- coreDat[nn.merge$value,tType]
    
    # get the real number of interaction
    realInt <- nrow(nn.merge.val[nn.merge.val$source >= 0.5 & nn.merge.val$value >= 0.5,])
    
    # permutation test, randomize the distribution of target 
    # get positions of target cell type
    
    # first, find out which cells are positive for source marker?
    sPos <- coreDat[which(coreDat[,sType] >= 0.5), c('Xcoord', 'Ycoord')]
    
    # second, find out which cells are positive for target marker?
    tPos <- coreDat[which(coreDat[,tType] >= 0.5), c('Xcoord', 'Ycoord')]
    
    
    # generate simulations of CSR
    simpp <- poisp(tPos, 50)
    
    #if(nrow(S))
    
    null <- lapply(simpp, function(x) {nnIntrxn(sPos, cbind(x$x, x$y))})
    
    nullDF <- data.frame(unlist(null))
    colnames(nullDF) <- 'simulates'
    
    
    #jpeg('./Figures/Simulations_35.jpeg', unit = 'in', width = 10, height = 5, res = 300)
    #ggplot(nullDF, aes(x = simulates)) + 
    #  theme_void() +
    #  geom_histogram(aes(y =..density..),
    #                 bins = 30,
    #                 colour = "black", 
    #                 fill = "#e4e1e6") +
    #  stat_function(fun = dnorm, args = list(mean = mean(nullDF$simulates), sd = sd(nullDF$simulates)), size = 1)
    #dev.off()
    #temp <- unlist(null)
    
    z.score <- (realInt - mean(unlist(null))) / sd(unlist(null))
    
    Heatmap[as.character(grid[id,1]), as.character(grid[id, 2])] <- z.score
    
    
  }
  

  print(core)
  write.csv(Heatmap, paste('./Data/Zscores/Z_score_', core ,'v2.csv', sep = ''))
  
  
}


#res.pca <- prcomp(textureFeature_all, scale = TRUE)

#res <- kmeans(textureFeature_all, centers = 3)

#res$cluster


#fviz_nbclust(textureFeature_all, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot")



plot(res.pca$x)

for(sType in mType){
  
  
  #print(rid)
  #sType <- 'CD8a'
  cid <- 1
  
  for(tType in mType){
    
    #print(cid)
    #tType <- 'CD4'
    # For each pair in nn.merge, find all rows that both elements are positive 
    # for this marker
    
    # expression for source column
    
    nn.merge.val <- nn.merge
    
    nn.merge.val$source <- coreDat[nn.merge$source,sType]
    nn.merge.val$value <- coreDat[nn.merge$value,tType]
    
    # get the real number of interaction
    realInt <- nrow(nn.merge.val[nn.merge.val$source >= 0.5 & nn.merge.val$value >= 0.5,])
    
    # permutation test, randomize the distribution of target 
    # get positions of target cell type
    
    # first, find out which cells are positive for source marker?
    sPos <- coreDat[which(coreDat[,sType] >= 0.5), c('Xcoord', 'Ycoord')]
    
    # second, find out which cells are positive for target marker?
    tPos <- coreDat[which(coreDat[,tType] >= 0.5), c('Xcoord', 'Ycoord')]
    
    
    # generate simulations of CSR
    simpp <- poisp(tPos, 500)
    
    #if(nrow(S))
    
    null <- lapply(simpp, function(x) {nnIntrxn(sPos, cbind(x$x, x$y))})
    
    nullDF <- data.frame(unlist(null))
    colnames(nullDF) <- 'simulates'
    
    
    #jpeg('./Figures/Simulations_35.jpeg', unit = 'in', width = 10, height = 5, res = 300)
    #ggplot(nullDF, aes(x = simulates)) + 
    #  theme_void() +
    #  geom_histogram(aes(y =..density..),
    #                 bins = 30,
    #                 colour = "black", 
    #                 fill = "#e4e1e6") +
    #  stat_function(fun = dnorm, args = list(mean = mean(nullDF$simulates), sd = sd(nullDF$simulates)), size = 1)
    #dev.off()
    #temp <- unlist(null)
    
    z.score <- (realInt - mean(unlist(null))) / sd(unlist(null))
    
    Heatmap[rid, cid] <- z.score
    
    
    cid <- cid + 1
    
  }
  
  
  
  