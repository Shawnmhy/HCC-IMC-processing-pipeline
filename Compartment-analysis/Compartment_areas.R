# Author: Haoyang Mi
# This script clean up the point patterns to find boundaries

############################
# Area:
# Core 6: 2735839622 + 2340245095 + 4283705590 + 383623784 + 2148503207 + 1731417348 (Tumor) -> 0.1945/mm2 (pixel width and height 264.58)
#         32025935612 - (2735839622 + 2340245095 + 4283705590 + 383623784 + 2148503207 + 1731417348) -> 0.263 / mm2
# Core 7: 26788000859 -> 0.383/mm2 (Tumor); 3210749066 -> 0.046 / mm2 (Immune)
# Core 8: 22448641818 -> 0.321/mm2 (Tumor); 6735467599 -> 0.096 / mm2 (Immune)
# Core 9: 27909330382 -> 0.399/mm2 (Tumor); 3727661115 -> 0.053 / mm2 (Immune)
#Core 22: 18536379258 -> 0.265/mm2 (Tumor); (30359832313 - 18536379258) -> 0.169 / mm2 (Immune)
#Core 23: 22785572708 -> 0.325/mm2 (Tumor); (29232692431 - 22785572708) -> 0.09 / mm2 (Immune)
#Core 27: (27159933919 - 1736247648) -> 0.363/mm2 (Tumor); 1736247648 -> 0.025 / mm2 (Immune)
#Core 28: 24323638069 -> 0.347/mm2 (Tumor); (27281111432 - 24323638069) -> 0.042 / mm2 (Immune)

(24323638069)/264.58/264.58/1000000
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

setwd("D:/DP/Projects/HCC/Compartments")

# write area to data

area <- data.frame(matrix(nrow = 8, ncol = 0))%>%
  mutate(core = c('6', '7', '8', '9', '22', '23', '27', '28'),
         areaImm.um = c(0.263, 0.046, 0.096, 0.053, 0.169, 0.092, 0.025, 0.042),
         areaHCC.um = c(0.195, 0.383, 0.321, 0.399, 0.265, 0.325, 0.363, 0.347))

write.csv(area, 'compartment_areas.csv', row.names = FALSE)
