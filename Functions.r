"
This document provide functions for spatial heterogenetiy analysis
@author: Haoyang Mi
"

library(spatstat)
library(dplyr)
library(flexclust)
library(plyr)
library(rjson)
library(spatstat)
library(pracma)



fac2num <- function(x){
  
  x <- as.numeric(as.character(x))
  return(x)
}

ShannonE <- function(types, coreData){

    #types <- ctype_tensor
    type.count <- length(types)

    Total <- nrow(coreDat) # total number of cells

    # all int/ext data
    ctype_stat_all <- data.frame(matrix(nrow =  0, ncol = 0))

    for (cseq in seq_len(type.count)) {   # outer loop, calcualte interior stats
        #cseq <- 1
        ctype_int <- 0
        ctype_ext <- 0
        p <- 0 # ratio

        # current cell type
        ctype <- types[cseq]

        # coordinates data for the current core, current cell type

        ctype.Dat <- coreDat[coreDat$ctype_no == ctype, c('Xcoord', 'Ycoord')]

        p <- nrow(ctype.Dat)/Total

        # interior score
        ctype_int <- mean(as.matrix(dist(ctype.Dat[c('Xcoord', 'Ycoord')])))


        # other cell types
        ctype_other <- types[-cseq]

        ctype_stat <- cbind(p, ctype_int)

        for (cseq_other in ctype_other) {

            ctype_other.Dat <- coreDat[coreDat$ctype_no == cseq_other, c('Xcoord', 'Ycoord')]

            # exterior score

            ctype_ext <- ctype_ext +  mean(dist2(ctype.Dat, ctype_other.Dat))

        }
        # number of computations = No. cell types - 1
        ctype_stat <- cbind(ctype_stat, ctype_ext / (type.count - 1))

        colnames(ctype_stat) <- c('p', 'int', 'ext')

        ctype_stat_all <- rbind(ctype_stat_all ,cbind(ctype_stat, ctype))

    }

    ShannonH <- 0
    # combine row data


    for (dat in seq_len(type.count)) {
      p <- as.numeric(as.character(ctype_stat_all[dat, 1]))
      d_int <- as.numeric(as.character(ctype_stat_all[dat, 2]))

      d_ext <- as.numeric(as.character(ctype_stat_all[dat, 3]))
      d_final <- d_int / d_ext

      if(isTRUE(d_int*d_ext == 0)){
        d_final <- 0
      }
      if(isTRUE(p != 0)){
        ShannonH <- -d_final*p*log2(p) + ShannonH
      }
    }

    return(ShannonH)
}



# convet graph to covariance matrix
#@ input: g -> graph 
#         null_matrix -> an empty matrix

graph2covMatrix <- function(Community_edges, node_type){
  
  
  from_types <- node_type[match(Community_edges$from, node_type$X1),]
  to_types <- node_type[match(Community_edges$to, node_type$X1),]
  
  nt_types <- data.frame(cbind(as.character(from_types$X2), as.character(to_types$X2)))
  
  # create covariance matrix
  covar_Matrix <- matrix(nrow = 17, ncol = 17)
  
  colnames(covar_Matrix) <- ctype_names
  rownames(covar_Matrix) <- ctype_names
  
  for(row in seq_len(17)){
    
    # from
    row_name <- ctype_names[row]
    for (col in seq_len(17)) {
      
      col_name <- ctype_names[col]
      covar_Matrix[row, col] <- nrow(nt_types[nt_types$X1 == row_name & nt_types$X2 == col_name,])+
        nrow(nt_types[nt_types$X1 == col_name & nt_types$X2 == row_name,])
      
    }
  }
  return(covar_Matrix)
}

GNN <- function(Community_edges, NodeFeature){

  
  
  # re-index
  NodeFeature$reindex <- seq_len(nrow(NodeFeature)) - 1
  edges <- list()
  features <- list()
  
  
  ## merge
  from_oldNode <- data.frame(Community_edges$from)
  colnames(from_oldNode) <- 'id'
  to_oldNode <- data.frame(Community_edges$to)
  colnames(to_oldNode) <- 'id'
  
  newNode <- data.frame(cbind(NodeFeature$members, NodeFeature$reindex))
  colnames(newNode) <- c('id', 'new_id')
  
  from_new <- vector()
  for(row in seq_len(nrow(from_oldNode))){
    
    old <- from_oldNode$id[row]
    new <- newNode[newNode$id == old, 2]
    
    from_new <- c(from_new, new)
  }
  
  to_new <- vector()
  for(row in seq_len(nrow(to_oldNode))){
    
    old <- to_oldNode$id[row]
    new <- newNode[newNode$id == old, 2]
    
    to_new <- c(to_new, new)
  }
  
  newEdge <- data.frame(cbind(from_new, to_new))
  # edges to list
  
  for(lid in seq_len(nrow(newEdge))) {
    edges[[lid]] <- as.numeric(newEdge[lid,])
  }
  
  F <- NodeFeature[,c(1,68, 69)]
  F$ctype_int <- as.character(F$ctype_int)
  
  # features to list
  
  for(fid in seq_len(nrow(NodeFeature))) {
    features[[fid]] <- F[fid,]$ctype_int 
    names(features)[fid] <- F[fid, 3]
  }
  
  # combine features and edges to one list
  list_final <- list(edges = edges, features = features)
  return(list_final)
}

cor.network <- function(gct_file, allCore){
  
  #allCore <- NR_core
  R_intercor <- data.frame(matrix(nrow = 0, ncol = 0))
  
  for (core in allCore) {
    
    #core <- 2
    # core - dendrogram
    Core_dendro <- gct_file[gct_file$id == core,]
    Core_dendro <- Core_dendro[order(Core_dendro$Community_id, decreasing = FALSE),]
    
    GNN_vec <- read.csv(paste('./Data/Features/', core, '/nci.csv', sep = ''), row.names = 1)
    intercor <- outer(1:nrow(GNN_vec),1:nrow(GNN_vec), FUN = Vectorize( function(i,j) cor.test(t(GNN_vec[i,]), t(GNN_vec[j,]), method = 'spearman')$p.value) )
    
    colnames(intercor) <- Core_dendro$dendrogram_cut
    rownames(intercor) <- Core_dendro$dendrogram_cut
    diag(intercor) <- -10
    
    # melt matrix
    melt_intercor <- melt(intercor)
    melt_intercor <- melt_intercor[!(melt_intercor$value == -10),]
    R_intercor <- rbind(R_intercor, melt_intercor)
    
  }
  
  Cor_network <- data.frame(matrix(nrow = 0, ncol = 0))
  for(row in seq_len(8)){
    #row <- 1
    #col <- 1
    for(col in row:8){
      
      dat <- R_intercor[R_intercor$Var1 == row & R_intercor$Var2 == col,]
      
      if(nrow(dat) >= 5){
        strength <- nrow(dat[dat$value < 0.05,])/nrow(dat)
        
      } else {
        strength <- 0
      }
      #hist(dat$value)
      
      #value <- p.adjust(dat$value, method = 'fdr')
      
      # normalized communication strength
      #strength <- length(value[value < 0.05])/length(value)
      Cor_network <- rbind(Cor_network, cbind(row, col, strength))
      
    }
  }
  Cor_network <- Cor_network[complete.cases(Cor_network$strength),]
  
  
  Cor_network <- Cor_network[Cor_network$strength > 0,]
 
  #test <- graph_from_data_frame(Cor_network[,1:2])
  #plot(test) 
 
  return(Cor_network)
}


#######################
# Kcross function #####
#######################



bivarAnalysis.Kcross <- function(type1, type2){
  
  #type1 <- typeA
  #type2 <- typeB
  
  colnames(type1) <- c('x', 'y')
  
  colnames(type2) <- c('x', 'y')
  
  if(nrow(type1)*nrow(type2) != 0){
    
    # read pts dat

    type1$attr <- 'ctypeA'
    
    type2$attr <- 'ctypeB'
    
    # create multitype df
    pts_OI <- rbind(type1, type2)
    
    # define the type
    species <- factor(pts_OI$attr)
    
    # create multitype ppp

    
    # check if empty  
    ppp1 <- ppp(type1$x, type1$y, owin(c(0,800), c(0, 800)))
    ppp2 <- ppp(type2$x, type2$y, owin(c(0,800), c(0, 800)))
    
    # prevent NA 
    
    if(is.empty(ppp1) == 'FALSE' & is.empty(ppp2) == 'FALSE'){
      
      
      
      
      multitype_ppp <- ppp(pts_OI$x, pts_OI$y, marks = species, owin(c(0, 800), c(0, 800)))
      K.cross <- data.frame(Kcross(multitype_ppp, i = 'ctypeA', j = 'ctypeB', r = seq(0,20,0.1), correction = 'Ripley'))
      
      #plot(Gihc)
      # relocat DF
      
      K.cross <- K.cross[complete.cases(K.cross),]
      
      
      
      # calculate the area (positive - negative )  
      K.cross$rs <- K.cross$iso - K.cross$theo
      
      A.to.B.diff.area <- trapz(K.cross$r, K.cross$iso) 
      
      
      # j to i
      
      K.cross <- data.frame(Kcross(multitype_ppp, i = 'ctypeB', j = 'ctypeA', r = seq(0,20,0.1), correction = 'Ripley'))
      
      K.cross <- K.cross[complete.cases(K.cross),]
      
      
      
      K.cross$iso <- K.cross$iso - K.cross$theo
      K.cross <- K.cross[complete.cases(K.cross),]
      B.to.A.diff.area <- trapz(K.cross$r, K.cross$iso) 
      
    } else {
      
      A.to.B.diff.area <- NA
      B.to.A.diff.area <- NA
      
    }
  } else {
    A.to.B.diff.area <- NA
    B.to.A.diff.area <- NA
  }
  return(list(A.to.B.diff.area, B.to.A.diff.area))
}



#----------------------------------------------------#
# generate poisson distribution ppp from a given ppp #
#----------------------------------------------------#


poisp <- function(pos_Target, nsim){
  
  #pos_Target <- tPos
  
  n <- nrow(pos_Target)
  
  simpp <- runifpoint(n = n, win = owin(c(0, 800), c(0, 800)), nsim = nsim)
  
  # return the result 
  return(simpp)
}


#-----------------------------#
# Find number of interactions #
#-----------------------------#

nnIntrxn <- function(sDF, simDF){ # sDF: the coordiantes for the source point pattern
  
  #sDF <- sPos
  if(nrow(sDF) == 0 | nrow(simDF) == 0){
    simInt <- 0
  } else{
    dist <- nn2(simDF, query = sDF, k = nrow(simDF), treetype = 'kd', searchtype = 'radius', radius = 20)
    
    nn.idx <- data.frame(dist$nn.idx)
    
    nn.idx$source <- seq_len(nrow(sDF))
    
    nn.merge <- melt(nn.idx, id.vars = 'source')
    
    # clear non-valid rows
    nn.merge <- nn.merge[nn.merge$value != 0, -2]
    
    
    # remove itself
    
    
    nn.merge$diff <- nn.merge$source - nn.merge$value
    nn.merge <- nn.merge[nn.merge$diff != 0, -3]
    
    simInt <- nrow(nn.merge)
  }

  
  # replace by real cell type
  return(simInt)
}



#-----------------------------#
#           cor.mtest         #
#-----------------------------#


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
}