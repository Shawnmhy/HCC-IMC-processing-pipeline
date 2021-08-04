# Network communications analysis #
# @Author: Haoyang Mi

library(ggplot2)
library(reshape2)
library(igraph)
library(reshape2)
library(tidyr)

setwd("D:/DP/Projects/HCC")

source('D:/DP/Projects/HCC/Functions.r')


# read gct file
gct_file <- data.frame(read.delim("D:/DP/Data/HCC/Community_Clustering.txt"))

# Bulk tumor community
gct_file$dendrogram_cut <- c(rep(1, 23), rep(2, 33), rep(3, 25), rep(4, 221), rep(5, 23), rep(6, 97), rep(7, 460), rep(8, 88))

# Responder
HCCdata <- readRDS("D:/DP/Data/HCC/hccdataset")

#remove UA Noncell

R_core <- unique(HCCdata[HCCdata$response == 'R',]$Core)
NR_core <- unique(HCCdata[HCCdata$response == 'NR',]$Core)


Rnet <- cor.network(gct_file, R_core)
NRnet <- cor.network(gct_file, NR_core)



#routes_network <- layout_components(routes_network)

g <- graph_from_data_frame(Rnet[,1:2])
plot(g)


# --------- Individual network count in R and NR -------------#

summary_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(core in seq_len(37)){
  #core <- 1
  summary <- gct_file %>%
    filter(id == core) %>%
    group_by(dendrogram_cut) %>%
    tally() %>%
    cbind(core)
  summary_all <- rbind(summary_all, summary)
  
  
  #print(unique(HCCdata[HCCdata$Core == core, 'response']))
  
}

# 
summary_all <- dcast(summary_all, core ~ dendrogram_cut, value.var = 'n')
summary_all[is.na(summary_all)] <- 0



# merge patient
Patient_table <- read.csv('Patient_Table.csv')
colnames(Patient_table)[1] <- 'core'
summary_all_patient <- merge(Patient_table, summary_all, by = 'core')


# R versus NR

R_count <- colSums(summary_all_patient[summary_all_patient$response == 'R', 4:11])
NR_count <- colSums(summary_all_patient[summary_all_patient$response == 'NR', 4:11])

response_count <- rbind(t(R_count), t(NR_count)) %>%
  data.frame()
row.names(response_count) <- c('R', 'NR')

require(tidyverse)
response_count <- response_count %>% rownames_to_column('group')

colnames(response_count) <- c('group', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
require(ggradar)

p <- ggradar(
  response_count[1:2, 1:8], 
  values.radar = c("0", "100", '300'),
  grid.min = 0, grid.mid = 100, grid.max = 300,
  group.line.width = 2, 
  group.point.size = 5,
  group.colours = c("#ef776d", "#21b7bd"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 7,
  grid.label.size = 8,
  )
p
ggsave(p, file=paste0("D:/DP/Projects/HCC/Figures/RadarPlot.png"), width = 8, height = 8, units = "in", dpi = 300)


write.csv(summary_all_patient, 'community_count_each_core.csv', row.names = FALSE)
