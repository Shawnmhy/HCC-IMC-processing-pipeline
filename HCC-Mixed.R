library(ggplot2)
library(reshape2)
library(tidyr)
library(matrixStats)
library(ggrepel)
library(gg.gap)
library(dplyr)
library(rstatix)
library(ggpubr)
library(descr)
library(pracma)
library(RANN)
library(flexclust)
library(ggvoronoi)

setwd("D:/DP/Projects/HCC")

source('D:/DP/Projects/HCC/Functions.r')

HCCdata <- readRDS("D:/DP/Data/HCC/hccdataset")
unique(HCCdata$ctype)

#HCCdata <- HCCdata%>%
#  filter(ctype != 'UA Noncell')

# rename Lag3 to LAG3

colnames(HCCdata)[11] <- 'LAG3'

# read core response data
res <- read.csv('D:/DP/Projects/HCC/Patient_Table.csv')
colnames(res)[1] <- 'Core'
# scale using quantile
rng <- colQuantiles(as.matrix(HCCdata[,1:32]), probs = c(0.01, 0.99))


HCCdata0 <- t((t(HCCdata[,1:32]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
HCCdata0[HCCdata0 < 0] <- 0; HCCdata0[HCCdata0 > 1] <- 1

HCCdata[,1:32] <- HCCdata0


# assign unique ID to each cell
HCCdata$uniqueID <- seq_len(nrow(HCCdata))


#-------------- Volcano plot ---------------#



exprCell <- HCCdata %>%
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  select(Core, GranB, LAG3, PDL1, Arg1, HLADR, CD28, CCR6, Ki67, Xcoord, Ycoord, ctype, uniqueID, response) %>%
  gather(key = 'marker', value = 'expr', -c('ctype', 'response', 'Core', 'Xcoord', 'Ycoord', 'uniqueID')) #%>%
  #group_by(ctype, marker) %>%
  #{
   # print(unique(with(response)))
    #if(length(unique(response)) > 1) 
    #wilcox_test(expr ~ response)
  
  #}

test_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(cid in unique(exprCell$ctype)){
  for(mkr in unique(exprCell$marker)){
    
    
    dat <- exprCell %>%
      filter(ctype == cid,
             marker == mkr)
    
    # NR count?
    NR_count <- nrow(dat[dat$response == 'NR',])
    R_count <- nrow(dat[dat$response == 'R',])
    
    
    if(NR_count > 100 & R_count > 100){
      
      test <- dat %>%
        wilcox_test(expr ~ response) %>%
        mutate(marker = mkr,
               ctype = cid)
      
      test$log2Fold <-  log2(mean(dat[dat$response == 'R', 'expr']) / mean(dat[dat$response == 'NR', 'expr']))
      test_all <- rbind(test_all, test)
    }
    
  }
}

test_all <- test_all %>%
  mutate(p.adjust = p.adjust(test_all$p, method = 'fdr')) %>%
  mutate(p.log = -log10(p.adjust)) %>%
  mutate(label = 1*(p.adjust < 0.01)) %>%
  mutate(annotation = paste0(marker, ' on ', ctype)) %>%
  filter(ctype != 'UA Noncell')

test_all$label <- as.factor(test_all$label)



# tailor markers for specific cell types

# add a global index
test_all$gid <- seq_len(nrow(test_all))

# all gid to exlude
id_to_exclude_all <- data.frame(matrix(nrow = 0, ncol = 0))
  
for(mkr in unique(test_all$marker)){
    
  # which row have this pattern? e.g. CCR6-, Arg1-
  idx <- grep(paste0(mkr, '-'), test_all$ctype, fixed = TRUE)
  
  id_to_exclude <- test_all[idx,] %>%
    filter(marker == mkr) %>%
    select(gid) %>%
    data.frame()
  
  id_to_exclude_all <- rbind(id_to_exclude_all, id_to_exclude)
}



# remove id_to_exclude_all from test_all

test_all <- test_all[-id_to_exclude_all$gid, ]


# select top 10 for each category

test_all_ranked <- test_all %>%
  group_by(label) %>% 
  mutate(rank = rank(p.adjust, ties.method = "first")) 
  
test_all_ranked[test_all_ranked$rank > 10 | test_all_ranked$label == '0', 'annotation'] <- NA



p <- ggplot(test_all_ranked, aes(log2Fold, p.log, label = annotation))+
  geom_point(shape = 21, aes(fill = label), size = 2) +
  theme_bw() +
  scale_fill_manual(values = c('0' = 'grey', '1' = 'red')) +
  xlab('Log2 Fold Difference (R/NR)') +
  ylab('-Log10 FDR-adjusted P-value') +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_blank()) +
  geom_text_repel(point.padding = unit(1, "lines"), color = ifelse(test_all_ranked$log2Fold < 0, '#d83637', '#0f3f80'))

#p <- gg.gap(p, segments = c(80, 100), ylim = c(0,300), tick_width = 20, c(0.8, 0, 0.2)) 

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/volcano.png"), width = 6, height = 7, units = "in", dpi = 300)


# -------------------------------- #
# Based on the previous volcano plot
# plot boxplot for furhter comparison
Marker <- 'Arg1'
Ctype <- 'Imm Mac (CD163-)'

mkrCell <- exprCell %>%
  filter(marker == Marker,
         ctype == Ctype)


stat.test <- mkrCell %>%
  wilcox_test(expr ~ response) %>%
  add_significance()
stat.test 

mkrCell$response <- factor(mkrCell$response, levels = c("R", "NR"))

p <- ggplot(mkrCell, aes(x = response, y = expr), fill = 'transparent') +
  geom_violin(position = 'dodge',  aes(color = response)) +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(mkrCell$expr) + 0.1, label.size = 15) +
  geom_boxplot(aes(ymin = min(expr),ymax = max(expr), color = response), fill = 'transparent') +
  theme(axis.title = element_text(size = 42),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 40),
        axis.title.x = element_blank(),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white', size = 2),
        #panel.spacing = unit(0, "mm"),
        axis.line = element_blank(),
        # remove spacing between facets
        #panel. = element_blank(),
        strip.text = element_text(size = 18),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  ylab("Normalized Expression" ) +
  scale_color_manual(values = c('R' = '#0f3f80', 'NR' = '#d83637'))
p
#ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/Arg1_CD163-Macrophage.png"), width = 6, height = 6, units = "in", dpi = 300)


# Histogram
Marker <- 'CCR6'
Ctype <- 'Imm Mac (CD163-)'
Thresh1 <- quantile(HCCdata[HCCdata$ctype == 'Imm Mac (CD163-)', 'Arg1'], probs = c(0.8))[1]
Thresh2 <- quantile(HCCdata[HCCdata$ctype == 'Imm Mac (CD163-)', 'CCR6'], probs = c(0.8))[1]

Density <- exprCell %>%
  filter(marker %in% c('CCR6', 'Arg1')) %>%
  filter(ctype == Ctype)


p <- ggplot(Density, aes(x=expr, color = marker)) +
  theme_bw() +
  geom_density(position="identity", size = 1) +
  scale_color_manual(values = c('Arg1' = '#99ccff', 'CCR6' = '#ffcc99')) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_blank(),
        legend.position = 'none') +
  geom_vline(aes(xintercept = Thresh1), linetype = 'dashed', size = 1, color = '#99ccff') +
  geom_vline(aes(xintercept = Thresh2), linetype = 'dashed', size = 1, color = '#ffcc99')

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/density.png"), width = 8, height = 6, units = "in", dpi = 300)

#-----------------------------------#
# May 3rd
# Based on signals discover above, can we use that to stratify cores?
#--------------------------#
#Marker <- 'CCR6'
#Ctype <- 'Imm Mac (CD163-)'
#--------------------------#
# first get the quantiles

for(num in seq(from = 0.5, to = 0.99, by = 0.01)){
  
  num <- 0.8
  CCR6_Thresh <- quantile(HCCdata[HCCdata$ctype == 'Imm Mac (CD163-)', 'CCR6'], probs = c(0.25, num))[2]
  
  Arg1_Thresh <- quantile(HCCdata[HCCdata$ctype == 'Imm Mac (CD163-)', 'Arg1'], probs = c(0.25, num))[2]
  
  
  M2Macro_ratio <- data.frame(matrix(nrow = 0, ncol = 0))
  scDF <- data.frame(matrix(nrow = 0, ncol = 0))
  for(core in unique(exprCell$Core)){
    
    # CD163- Macrophage count
    M2Macro <- nrow(exprCell[exprCell$ctype == 'Imm Mac (CD163-)' & exprCell$Core == core,])
    
    # response
    r <- unique(exprCell[exprCell$Core == core, 'response'])
    
    
    # get coordinates etc for selected cells
    corescDF <- HCCdata %>%
      filter(Core == core) %>% #Imm low: 19, 20, 21, 35, 37
      filter(ctype == 'Imm Mac (CD163-)') %>%
      select(Core, Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response) %>%
      filter(Arg1 > Arg1_Thresh | CCR6 > CCR6_Thresh)
    
    # selected count
    sc <- nrow(corescDF)
    
    
    M2Macro_ratio <- rbind(M2Macro_ratio, cbind(core, sc, sc/M2Macro, r))
    
    scDF <- rbind(scDF, corescDF)
  }
  
  M2Macro_ratio$sc <- as.numeric(as.character(M2Macro_ratio$sc))
  tst <- wilcox.test(M2Macro_ratio[M2Macro_ratio$r == 'NR', 'sc'], M2Macro_ratio[M2Macro_ratio$r == 'R', 'sc'])
  print(paste0('threshold: ', num,' pval = ', tst$p.value))
}



stat.test <- M2Macro_ratio %>%
  wilcox_test(sc ~ r, c('NR', 'R')) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "r")

M2Macro_ratio$r <- factor(M2Macro_ratio$r, levels = c("R", "NR"))


p <- ggboxplot(M2Macro_ratio, x = "r", y = "sc", add.params = list(size = 4),
               color = 'r', add = 'jitter', palette = c("#0f3f80", "#d83637"),
               bxp.errorbar = TRUE, size = 1) +
  stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(M2Macro_ratio$sc) + 10, label.size = 14) +
  
  theme(axis.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none', 
        axis.text = element_text(size = 24))

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/Hazard_count_RNR.png"), width = 5, height = 6, units = "in", dpi = 300)


#----------------------------#
# May 3rd
# Correlation plot for CCR6 + Arg1

cordata <- HCCdata %>%
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  filter(ctype == 'Imm Mac (CD163-)') %>%
  filter(response == 'NR') %>%
  select(Arg1, CCR6) %>%
  arrange(Arg1)

cor(cordata$Arg1, cordata$CCR6)
lm_model <- lm(CCR6 ~ Arg1, data = cordata)
lm_model
summary(lm_model)$r.squared

p <- ggplot(data = cordata, aes(Arg1, CCR6)) +
  theme_bw() +
  geom_point(shape = 21, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  coord_equal(ratio = 1) +
  xlim(0,0.9) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white', size = 2),
        #panel.spacing = unit(0, "mm"),
        axis.line = element_blank(),
        # remove spacing between facets
        #panel. = element_blank(),
        strip.text = element_text(size = 18),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/Cor-Arg1-CCR6-NR.png"), width = 6, height = 6, units = "in", dpi = 300)



cordata <- HCCdata %>%
  filter(ctype == 'Imm Mac (CD163-)') %>%
  select(Arg1, CCR6) %>%
  arrange(Arg1)


p <- ggplot(data = cordata, aes(Arg1, CCR6)) +
  theme_bw() +
  geom_point(shape = 21, size = 1, color = 'grey') +
  geom_smooth(method = "lm", se = TRUE) +
  coord_equal(ratio = 1) +
  xlim(0,1) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'white', size = 2),
        #panel.spacing = unit(0, "mm"),
        axis.line = element_blank(),
        # remove spacing between facets
        #panel. = element_blank(),
        strip.text = element_text(size = 18),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/Cor-Arg1-CCR6.png"), width = 6, height = 6, units = "in", dpi = 300)



#---------------------------------------#
# May 4th
# neighborhood analysis of hazard cells

# remove self points
HCCdata_subset <- HCCdata %>%
  filter(!(uniqueID %in% scDF$uniqueID))

scNN_all <- data.frame(matrix(nrow = 0, ncol = 0))
non_scNN_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(core in unique(exprCell$Core)){
  
  scDF_core <- scDF %>%
    filter(Core == core)
  
  HCCdata_subset <- HCCdata %>%
    filter(!(uniqueID %in% scDF$uniqueID)) %>%
    filter(Core == core)
  
  # Hazard CD163- Macrophages neighboring cells
  
  scNN <- nn2(data = HCCdata_subset[, c('Xcoord', 'Ycoord')], query = scDF_core[, c('Xcoord', 'Ycoord')], k = nrow(HCCdata_subset),
      treetype = 'kd', searchtype = 'radius', radius = 20) %>%
    with(nn.idx) %>%
    data.frame() %>%
    #select(-c('X1')) %>% # remove identity
    mutate(total = rowSums(across()), id = seq_len(nrow(scDF_core))) %>% # row sum to find no neighbor cells
    filter(total > 0) %>% # remove no neighbors cells
    select(-'total') %>%
    gather(key = 'source', value = 'target', -'id') %>%
    filter(target > 0) %>% # remove no neighbors cells
    select(-'source') %>%
    `colnames<-`(c("source", "target")) %>%
    mutate(M2Mac = pmin(source, target), # remove duplicates
           nn = pmax(source, target)) %>%
    distinct(M2Mac, nn)
  
  # Non-Hazard CD163- Macrophages neighboring cells
  non_scDF <- HCCdata %>%
    filter(!(uniqueID %in% scDF$uniqueID)) %>%
    filter(ctype == 'Imm Mac (CD163-)') %>%
    filter(Core == core) %>%
    select(Core, Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response)
  
  HCCdata_subset2 <- HCCdata_subset %>%
    filter(!(uniqueID %in% non_scDF$uniqueID))
  
  
  non_scNN <- nn2(data = HCCdata_subset2[, c('Xcoord', 'Ycoord')], query = non_scDF[, c('Xcoord', 'Ycoord')], k = nrow(HCCdata_subset2),
              treetype = 'kd', searchtype = 'radius', radius = 20) %>%
    with(nn.idx) %>%
    data.frame() %>%
    #select(-c('X1')) %>% # remove identity
    mutate(total = rowSums(across()), id = seq_len(nrow(non_scDF))) %>% # row sum to find no neighbor cells
    filter(total > 0) %>% # remove no neighbors cells
    select(-'total') %>%
    gather(key = 'source', value = 'target', -'id') %>%
    filter(target > 0) %>% # remove no neighbors cells
    select(-'source') %>%
    `colnames<-`(c("source", "target")) %>%
    mutate(M1Mac_neg = pmin(source, target), # remove duplicates
           nn = pmax(source, target)) %>%
    distinct(M1Mac_neg, nn)
  
  
  
  
  
  
  
  scNN$source <- 'Imm Mac (CD163-)_Pos'
  scNN$target <- HCCdata[scNN$nn, 'ctype']
  scNN$response <- unique(scDF_core$response)
  scNN$Core <- core
  
  
  
  non_scNN$source <- 'Imm Mac (CD163-)_Neg'
  non_scNN$target <- HCCdata[non_scNN$nn, 'ctype']
  non_scNN$response <- unique(scDF_core$response)
  non_scNN$Core <- core
  
  
  
  
  scNN_all <- rbind(scNN_all, scNN)
  non_scNN_all <- rbind(non_scNN_all, non_scNN)
  
}

# each cell type counts in R


non_scNN_tb <- data.frame(table(non_scNN_all[, 'target'])) %>%
  filter(Var1 != 'UA Noncell')
scNN_tb <- data.frame(table(scNN_all[, 'target'])) %>%
  filter(Var1 != 'UA Noncell')




scNN_paried <- merge(non_scNN_tb, scNN_tb, by = 'Var1') %>%
  `colnames<-`(c('ctype', 'Neg', 'Pos')) %>%
  filter(ctype != 'Imm Mac (CD163-)') %>%
  mutate(ratio_Neg = Neg / sum(Neg),
         ratio_Pos = Pos / sum(Pos),
         diff = ratio_Pos - ratio_Neg)


ggpaired(scNN_paried, cond1 = "ratio_Pos", cond2 = "ratio_Neg",
         fill = "condition", palette = "jco")






#---------------------------------------#
# May 4th
# Distance analysis of hazard cells



PO_dist_all <- data.frame(matrix(nrow = 0, ncol = 0))
NO_dist_all <- data.frame(matrix(nrow = 0, ncol = 0))

# Distance from positive Macrophages to other cell types

for(core in unique(exprCell$Core)){
  
  # All Cell Types Other than macrophages
  
  Other <- HCCdata %>%
    filter(Core == core) %>%
    filter(ctype != 'Imm Mac (CD163-)') %>%
    select(Xcoord, Ycoord, ctype, uniqueID, response)
  
  
  
  # All 'Positive' Macrophages (CD163-)
  
  posMac <- HCCdata %>%
    filter(Core == core) %>%
    filter(ctype == 'Imm Mac (CD163-)') %>%
    select(Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response) %>%
    filter(Arg1 > Arg1_Thresh | CCR6 > CCR6_Thresh)
  
  
  # All 'Negative' Macrophages (CD163-)
  
  negMac <- HCCdata %>%
    filter(Core == core) %>%
    filter(ctype == 'Imm Mac (CD163-)') %>%
    select(Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response) %>%
    filter(!(uniqueID %in% posMac$uniqueID))
  
  
  # Distance between Pos to Other
  PO_dist <- dist2(Other[, c('Xcoord', 'Ycoord')], posMac[, c('Xcoord', 'Ycoord')]) %>%
    data.frame() %>%
    mutate(ctype = Other$ctype) %>% # add cell type
    gather(key = 'id', value = 'dist', -'ctype') %>% # gather
    select(ctype, dist)# remove X1, X2, ....
    
    
    # Distance between Neg to Other
  NO_dist <- dist2(Other[, c('Xcoord', 'Ycoord')], negMac[, c('Xcoord', 'Ycoord')]) %>%
    data.frame() %>%
    mutate(ctype = Other$ctype) %>% # add cell type
    gather(key = 'id', value = 'dist', -'ctype') %>% # gather
    select(ctype, dist)# remove X1, X2, ....
  
  
  
  PO_dist_all <- rbind(PO_dist_all, PO_dist)
  NO_dist_all <- rbind(NO_dist_all, NO_dist)
  
}


PO_dist <- PO_dist_all %>%
  group_by(ctype) %>%
  dplyr::summarize(meanDist = mean(dist))

NO_dist <- NO_dist_all %>%
  group_by(ctype) %>%
  dplyr::summarize(meanDist = mean(dist))


NOPO_paried <- merge(PO_dist, NO_dist, by = 'ctype') %>%
  `colnames<-`(c('ctype', 'PO_dist', 'NO_dist')) %>%
  mutate(`PO - NO` = (PO_dist - NO_dist) / NO_dist) %>%
  filter(ctype != 'UA Noncell')
  



stat.test <- NOPO_paried %>%
  gather(key = 'id', value = 'dist', -'ctype') %>%
  wilcox_test(dist ~ id, comparisons = list(c('PO_dist', 'NO_dist')), paired = TRUE) %>%
  add_significance()

#stat.test <- stat.test %>% add_xy_position(x = "id")



p <- ggpaired(NOPO_paried, cond1 = "NO_dist", cond2 = "PO_dist",
             color = "condition",  width = 0.5,
             palette = c("#74a4c7", "#e06868"),  point.size = 3,
             scale_x_discrete(labels = c("draft", "final")),
             legend = 'none', 
             ylab = expression(paste('Distance, ', mu, 'm')),
             xlab = FALSE, font.title = list(size = 20)) +
  font("xy.text", size = 18) +
  font('y.title', size = 20) +
  scale_x_discrete(labels= c(expression(paste('Non-Hazard M', Phi)), expression(paste('Hazard M', Phi))))
p  

ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/Hazard_Macro_compare.png"), width = 5, height = 5, units = "in", dpi = 300)



# all distances, not just shortest distance

posMac <- HCCdata %>%
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  filter(ctype == 'Imm Mac (CD163-)') %>%
  select(Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response) %>%
  filter(Arg1 > Arg1_Thresh | CCR6 > CCR6_Thresh)

# All 'Negative' Macrophages (CD163-)

negMac <- HCCdata %>%
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  filter(ctype == 'Imm Mac (CD163-)') %>%
  select(Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response) %>%
  filter(!(uniqueID %in% posMac$uniqueID))



Bcells <- HCCdata %>%
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  filter(ctype == 'Imm B') %>%
  select(Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response)




t <- scDF %>%
  group_by(Core) %>%
  tally()

#---------------------------------------#
# May 5th
# Voronoi tesselation: B cell and Hazard/Non-hazard macrophages

non_scDF <- HCCdata %>%
  filter(!(uniqueID %in% scDF$uniqueID)) %>%
  filter(ctype == 'Imm Mac (CD163-)') %>%
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  select(Core, Arg1, CCR6, Xcoord, Ycoord, ctype, uniqueID, response)


for(core in c(1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 25, 26, 29, 30, 31, 32, 33, 34, 36)){
  
  plotData <- HCCdata %>%
    filter(Core == core)
  
  Xmean <- (max(plotData$Xcoord) + min(plotData$Xcoord))/2
  Ymean <- (max(plotData$Ycoord) + min(plotData$Ycoord))/2
  #centroid(cbind(core$Xcoord, core$Ycoord))
  
  n <- 1000 # number of points you want on the unit circle
  pts.circle <- data.frame(t(sapply(1:n,function(r)c(380*cos(2*r*pi/n),380*sin(2*r*pi/n)))))
  
  pts.circle[,1] <- pts.circle[,1] + Xmean
  pts.circle[,2] <- pts.circle[,2] + Ymean
  
  
  plotData$altLabel <- '0'
  
  plotData[plotData$ctype == 'Imm CD8 T', "altLabel"] <- '1'
  
  plotData[plotData$uniqueID %in% scDF$uniqueID, "altLabel"] <- '2'
  
  #plotData[plotData$uniqueID %in% non_scDF$uniqueID, "altLabel"] <- '3'
  
  plotData[plotData$ctype == 'Imm Treg', "altLabel"] <- '3'
  
  plotData[plotData$ctype == 'Imm DP T', "altLabel"] <- '4'
  
  plotData[plotData$ctype == 'Imm Neut (PDL1+)', "altLabel"] <- '5'
  
  
  require(ggvoronoi)
  
  p <- ggplot(plotData,aes(Xcoord, Ycoord)) +
    theme_void() +
    theme(legend.position = 'none') +
    #xlim(200, 500) +
    # ylim(250, 500) +
    geom_voronoi(aes(fill = altLabel), color = 'black', outline = pts.circle) +
    scale_fill_manual(values = c('0' = '#f6f6f6', '1' = '#74a4c7', '2' = '#e06868', '3' = '#63a078',
                                 '4' = '#eae551', '5' = '#c26715'))  +
    #geom_text(data = core, aes(x = Xcoord, y = Ycoord, label = ctype_integer), size = 6) +
    coord_fixed(ratio = 1)
  p
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/", core, ".png"), width = 8, height = 8, units = "in", dpi = 300)
  
}





# SpatialScore

SpatialScore_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(core in unique(exprCell$Core)){
  
  
  Response <- unique(exprCell[exprCell$Core == core, 'response'])
  # CD4 T cell
  
  CD8T <- HCCdata %>%
    filter(Core == core) %>%
    filter(ctype == 'Imm CD8 T') %>%
    select(Xcoord, Ycoord, ctype, uniqueID, response)
  
  # Tumor cells
  Tumor <- HCCdata %>%
    filter(Core == core) %>%
    filter(ctype %in% c('Imm CD4 T')) %>%
    select(Xcoord, Ycoord, ctype, uniqueID, response)
  
  # CD8 cells
  HazardMac <- HCCdata %>%
    filter(Core == core) %>%
    filter(ctype == 'Imm Mac (CD163-)') %>%
    filter(Arg1 > Arg1_Thresh) %>%
    select(Xcoord, Ycoord, ctype, uniqueID, response)
  
  # Nearest Tumor Distance
  
  print(nrow(CD8T) * nrow(HazardMac) * nrow(Tumor))
  if(nrow(CD8T) > 0 & nrow(HazardMac) >0 & nrow(Tumor) > 0){
    nnT <- dist2(CD8T[, c('Xcoord', 'Ycoord')], Tumor[, c('Xcoord', 'Ycoord')]) %>%
      rowMins()
    
    # Nearest CD8 Distance
    
    nnH <- dist2(CD8T[, c('Xcoord', 'Ycoord')], HazardMac[, c('Xcoord', 'Ycoord')]) %>%
      rowMins()
    
    SpatialScore <- nnT/ (nnH + nnT)
    
    SpatialScore_all <- rbind(SpatialScore_all, cbind(SpatialScore, core, Response, CD8T$uniqueID))
    
  }

}




SpatialScore_all$SpatialScore <- as.numeric(as.character(SpatialScore_all$SpatialScore))
colnames(SpatialScore_all)[4] <- 'uniqueID'


# Density: Spatial Score

p <- ggplot(SpatialScore_all, aes(x= SpatialScore, color = Response)) +
  theme_bw() +
  geom_density(position="identity", size = 1) +
  scale_color_manual(values = c('R' = '#21b7bc', 'NR' = '#ef776d')) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_blank(),
        legend.position = 'none') +
  ylim(0, 2.5)

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/density_riskscore_split.png"), width = 8, height = 6, units = "in", dpi = 300)

p <- ggplot(SpatialScore_all, aes(x= SpatialScore)) +
  theme_bw() +
  geom_density(position="identity", size = 1) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_blank(),
        legend.position = 'none') +
  geom_vline(aes(xintercept = 0.3), linetype = 'dashed', size = 1, color = '#375579') +
  geom_vline(aes(xintercept = 0.7), linetype = 'dashed', size = 1, color = '#7f9ec0') +
  ylim(0, 2.5)

p

ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/density_riskscore.png"), width = 8, height = 6, units = "in", dpi = 300)





stat.test <- SpatialScore_all %>%
  wilcox_test(SpatialScore ~ Response, c('R', 'NR')) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Response")

SpatialScore_all$Response <- factor(SpatialScore_all$Response, levels = c("R", "NR"))


p <- ggboxplot(SpatialScore_all, x = "Response", y = "SpatialScore", add.params = list(size = 4),
               color = 'Response', palette = c("#0f3f80", "#d83637"), 
               bxp.errorbar = TRUE, size = 1) +
  stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(SpatialScore_all$SpatialScore) + 0.1, label.size = 10) +
  
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none', 
        axis.text = element_text(size = 24), 
        axis.title.y = element_text(size = 26,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  ylab(expression(paste('CD8'^"+", ' T Cell Risk Score')))


p

ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/SpatialScore.png"), width = 5, height = 6, units = "in", dpi = 300)


# Density plot

p <- ggplot(SpatialScore_all, aes(SpatialScore, fill = Response, colour = Response)) +
  geom_density(alpha = 0.1)



p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/SpatialScore.png"), width = 5, height = 6, units = "in", dpi = 300)


# Granzyme B expression on high spatial score CD8 T cells

#thresh <- quantile(SpatialScore_all$SpatialScore, probs = c(0.4, 0.6))
#thresh

expr_SpatialScore <- HCCdata %>%
  #filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  merge(SpatialScore_all, by = 'uniqueID') %>%
  mutate(Group = 1 * (SpatialScore >= 0.3) + 1 * (SpatialScore >= 0.7)) %>%
  select(GranB, Group)


# GranB density plot 
expr_SpatialScore$Group <- as.character(expr_SpatialScore$Group)
p <- ggplot(expr_SpatialScore, aes(x= GranB, color = Group)) +
  geom_density(position="identity", size = 1) +
  theme_bw() +
  scale_color_manual(values = c('0' = '#375579', '1' = '#8d959c', '2' = '#7f9ec0')) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_blank(),
        legend.position = 'none'
        ) +
  ylim(0, 2)

p

ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/density_CD28.png"), width = 8, height = 6, units = "in", dpi = 300)




stat.test <- expr_SpatialScore %>%
  wilcox_test(GranB ~ Group, comparisons = list(c('0', '1'), c('1', '2'))) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Group")

expr_SpatialScore$Group <- factor(expr_SpatialScore$Group, levels = c("0", "1", '2'))


p <- ggboxplot(expr_SpatialScore, x = "Group", y = "GranB", add.params = list(size = 4),
               color = 'Group', palette = c("#375579", "#8d959c", '#7f9ec0'), bxp.errorbar = TRUE, size = 1) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = max(expr_SpatialScore$GranB) + 0.05, label.size = 10) +
  
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none', 
        axis.text = element_text(size = 24), 
        axis.title.y = element_text(size = 26,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  ylab('Normalized Expressions')

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/Mixed/GranB_SpatialScore.png"), width = 6, height = 6, units = "in", dpi = 300)





# CD8+ T cell

Tc_CD4 <- HCCdata %>%
  filter(ctype %in% c('Imm CD4 T')) %>% # CD8 T cells
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  select(GranB, Lag3, CD8a, Core, response, Xcoord, Ycoord, uniqueID) %>%
  filter((Lag3 + GranB) !=0 ) %>%
  mutate(EffScore = GranB / (GranB + Lag3), group = 'CD4')
  #mutate(exhaustInd = 1 * (GranB < 0.205 & Lag3 > 0.571))

Tc_CD8 <- HCCdata %>%
  filter(ctype %in% c('Imm CD8 T')) %>% # CD8 T cells
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  select(GranB, Lag3, CD8a, Core, response, Xcoord, Ycoord, uniqueID) %>%
  filter((Lag3 + GranB) !=0 ) %>%
  mutate(EffScore = GranB / (GranB + Lag3), group = 'CD8')


Tc <- rbind(Tc_CD4, Tc_CD8)
# Rename Lag3 to LAG3
colnames(Tc)[2] <- 'LAG3'


# Test and Boxplot

#------ BLOCK 1: GranB on CD4 + CD8
{
  stat.test <- Tc %>%
    group_by(group) %>%
    wilcox_test(GranB ~ response, comparisons = list(c('NR', 'R'), c('NR', 'R'))) %>%
    add_significance()
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  Tc$response <- factor(Tc$response, levels = c("R", "NR"))
  
  p <- ggboxplot(Tc, x = "group", y = "GranB",
                 color = 'response', palette = c("#0f3f80", "#d83637"),  size = 1) +
    stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(Tc$GranB) + 0.05, label.size = 14) +
    
    theme(axis.title = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          #legend.position = 'none', 
          axis.text = element_text(size = 24))
  
  p
  
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/GranB_mixed_compare.png"), width = 8, height = 6, units = "in", dpi = 300)
  
}

#------ BLOCK 2: LAG3 on CD4 + CD8
{
  stat.test <- Tc %>%
    group_by(group) %>%
    wilcox_test(LAG3 ~ response, comparisons = list(c('NR', 'R'), c('NR', 'R'))) %>%
    add_significance()
  
  stat.test <- stat.test %>% add_xy_position(x = "group")

  
  Tc$response <- factor(Tc$response, levels = c("R", "NR"))
  
  p <- ggboxplot(Tc, x = "group", y = "LAG3",
                 color = 'response', palette = c("#0f3f80", "#d83637"),  size = 1) +
    stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(Tc$GranB) + 0.05, label.size = 10) +
    
    theme(axis.title = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          #legend.position = 'none', 
          axis.text = element_text(size = 24))
  
  p
  
  ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/LAG3_mixed_compare.png"), width = 8, height = 6, units = "in", dpi = 300)
  
}

# Set up threshold
quantile <- HCCdata %>%
  filter(ctype %in% c('Imm CD4 T')) %>% # CD8 T cells
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  select(GranB, Lag3) %>%
  as.matrix() %>%
  colQuantiles(probs = c(0.25, 0.75))

quantile
# by core
Tc_core <- Tc_CD4 %>%
  mutate(Score = 1 * (GranB > 0.429) + 1* (Lag3 < 0.177)) %>%
  group_by(Core, response) %>%
  dplyr::summarize(mean = mean(Score)) %>%
  as.data.frame()
  
wilcox.test(Tc_core[Tc_core$response == 'R', 'mean'], Tc_core[Tc_core$response == 'NR', 'mean'])


# Mean Score


stat.test <- Tc_core %>%
  wilcox_test(mean ~ response, comparisons = list(c('R', 'NR'))) %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "response")


p <- ggboxplot(Tc_core, x = "response", y = "mean", add = 'jitter', add.params = list(size = 3),
               color = 'response', palette = c("#0f3f80", "#d83637"),  size = 1) +
  stat_pvalue_manual(stat.test, label = "p.signif", y.position = max(Tc$GranB) + 0.05, label.size = 10) +
  
  theme(axis.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none', 
        axis.text = element_text(size = 24))

p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/EffScore.png"), width = 4, height = 6, units = "in", dpi = 300)




#---------------------#
# Correlation Cancer Cell counts and Score


CancerCell <- HCCdata %>%
  filter(ctype %in% c('Hep', 'Hep Prolif', 'Hep Apop (CCR6-Arg1-)', 'Hep Apop', 'Hep (CCR6-Arg1-)')) %>% # CD8 T cells
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  group_by(Core) %>% 
  tally()



CancerCell_Score <- merge(CancerCell, Tc_core, by = 'Core')



LAG3 <- HCCdata %>%
  filter(ctype %in% c('Imm CD4 T')) %>% # CD8 T cells
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>%
  select(Lag3, GranB)


cor(LAG3)

Tc_count <- HCCdata %>%
  filter(ctype %in% c('Imm CD8 T', 'Imm CD4 T')) %>% # CD8 T cells
  filter(!(Core %in% c(6, 7, 8, 9, 22, 23, 27, 28, 19, 20, 21, 35, 37))) %>% #Imm low: 19, 20, 21, 35, 37
  group_by(Core) %>%
  tally() %>%
  filter(Core %in% Tc_core$Core)


Tc_core <- Tc_core %>%
  mutate(proportion = n / Tc_count$n) %>%
  merge(res, by = 'Core')


# neighborhood analysis
for(core in unique(Tc$Core)){
  
  # Tc data
  exTc_subset <- Tc %>%
    filter(Core == core) %>%
    filter(exhaustInd == 1)
  
  # original data
  HCCdata_subset <- HCCdata %>%
    filter(Core == core) %>%
    filter(!(uniqueID %in% exTc_subset$uniqueID)) %>% # remove exhausted CD8 T cells 
    select(Xcoord, Ycoord, ctype_no, ctype)

  neighbors <- nn2(data = HCCdata_subset[,c('Xcoord', 'Ycoord')], query = exTc_subset[,c('Xcoord', 'Ycoord')], k = nrow(HCCdata_subset), treetype = 'kd', searchtype = 'radius', radius = 20) %>%
    with(nn.idx) %>%
    data.frame() %>%
    #select(-c('X1')) %>% # remove identity
    mutate(total = rowSums(across()), id = seq_len(nrow(exTc_subset))) %>% # row sum to find no neighbor cells
    filter(total > 0) %>% # remove no neighbors cells
    select(-'total') %>%
    gather(key = 'source', value = 'target', -'id') %>%
    filter(target > 0) %>% # remove no neighbors cells
    select(-'source') %>%
    `colnames<-`(c("source", "target")) %>%
    mutate(exTc = pmin(source, target), # remove duplicates
           nn = pmax(source, target)) %>%
    distinct(exTc, nn)
  
  table(HCCdata_subset[neighbors$nn, 'ctype'])

}


nrow(Tc[Tc$exhaustInd == 1 & Tc$response == 'R',])/nrow(Tc)
nrow(Tc[Tc$exhaustInd == 1 & Tc$response == 'NR',])/nrow(Tc)

plot(Tc[Tc$exhaustInd == 1 & Tc$response == 'R',])