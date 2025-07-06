# Function Documentation
# Optimized functions
########################
# Packages Required
########################

library(INLA) # For Imputation
library(gamlss.spatial) # For SNR
library(tidyverse) # Data manipulation
library(igraph)
library(doParallel)
library(scales)
library(genie)
library(ggraph)
library(ar.matrix)
library(umap)
library(progress)
library(viridis)

SpatialDENet <- function(all_data,
                         metadata,
                         Sample_id = "Sample_id",
                         joint=TRUE,
                         countmodel = "nbinomial",
                         CollectPostD = FALSE,
                         offset = c(1, 1),
                         max.edge = c(.7, .7),
                         cutoff =1){
  
 
  n_genes = nrow(all_data)
  un =  unique(metadata[,Sample_id])
  Metadata_sub = NULL
  All_data = NULL
  for (k in un) {
    
    metadata_sub = metadata[metadata$Sample_id==k,]
    All_data = cbind(All_data,all_data[,metadata$Sample_id==k])
    
    coord = metadata_sub[,c("x","y")]
    coord[,"x"] = as.numeric(scale(coord[,"x"]))
    coord[,"y"] = as.numeric(scale(coord[,"y"]))
    
    mesh = inla.mesh.2d(loc.domain = coord, offset = offset, 
                        max.edge = max.edge, cutoff =cutoff)
    dmesh = book.mesh.dual(mesh=mesh)
    dmesh2 = st_as_sf(dmesh)
    coordsPoints = st_as_sf(coord, coords =c(1,2))
    aux = st_intersects(coordsPoints, dmesh2)
    aux = aux %>% unlist()
    
    metadata_sub <- bind_cols(metadata_sub,meshid=paste0(aux,k))
    Metadata_sub = rbind(Metadata_sub,metadata_sub)
    
  }
  all_data = as.matrix(All_data)
  Metadata_sub$meshid_id = as.numeric(as.factor(Metadata_sub$meshid))
  Centers = bind_cols(col=Metadata_sub$meshid_id,t(all_data)) %>% group_by(col) %>% summarise_all(mean)%>%dplyr::select(-col)
  Res = data.frame(col=Metadata_sub$meshid_id) %>% group_by(col) %>%summarise(nn = n())
  
  mst = ClusterToTree(Centers)
  
  # Get network using only spatial coordinates
  ############################################
  
  # Extract edges from mesh
  #edges <- mesh$graph$tv  # triangle vertices
  
  # Number of nodes
  #n_nodes <- mesh$n
  
  # Initialize sparse adjacency matrix
  #adj <- Matrix(0, nrow = n_nodes, ncol = n_nodes, sparse = TRUE)
  
  # Loop through triangles and set adjacency
  #for (i in 1:nrow(edges)) {
  #  tri <- edges[i, ]
  #  adj[tri[1], tri[2]] <- 1
  #  adj[tri[2], tri[1]] <- 1
  #  adj[tri[1], tri[3]] <- 1
  #  adj[tri[3], tri[1]] <- 1
  #  adj[tri[2], tri[3]] <- 1
  #  adj[tri[3], tri[2]] <- 1
  #}
  
  #mst <- graph_from_adjacency_matrix(adj[sort(unique(aux)),sort(unique(aux))], mode = "undirected", diag = FALSE)
  
  
  #####################
  # create mesh
  ###################
  
  A <- as_adjacency_matrix(mst, sparse = TRUE)
  
  # 3. Compute degree matrix
  D <- Diagonal(x = rowSums(A))
  
  # 4. Compute precision matrix Q = D - A
  Q <- D - A
  Q_sparse <- as(Q, "dgTMatrix")
  n= dim(Q_sparse)[1]
  # Indices
  
  idx1 <- as.numeric(V(mst))
  idx2 <- idx1 + max(V(mst))  # Shifted index to differentiate between the two stacks
  
  Graphid = data.frame(idx1,idx2=idx1, meshid_id = sort(unique(Metadata_sub$meshid_id)))
  Metadata_sub = left_join(Metadata_sub,Graphid,by="meshid_id")
  #metadata$idx2[metadata$Group=="Treatment"] = metadata$idx2[metadata$Group=="Treatment"]+vcount(mst)
  
  ##### Model gene-by-gene
  
  pvalue = sig_cont = sig_tactv = SNR= rep(0, n_genes)
  names(pvalue) = names(sig_cont) = names(sig_tactv) = names(SNR) = rownames(all_data)
  pvalue = SNR = matrix(NA,2,n_genes)
  colnames(pvalue) =colnames(SNR)  = rownames(all_data)
  treeEffect_actv = treeEffect_cont = matrix(NA, nrow =vcount(mst),ncol = n_genes)
  colnames(treeEffect_actv) = colnames(treeEffect_cont) = rownames(all_data)
  PosteriorDist = list()
  
  pb <- progress_bar$new(total = n_genes)
  
  for (j in 1:n_genes) {
    pb$tick()
    if(joint){
    y1 = ifelse(all_data[j,]>0,all_data[j,],NA)
    y2 = ifelse(all_data[j,]==0,0,1)
    n = length(y1)
    # Response list: each element of the list is a vector of length n * 2 (for 2 likelihoods)
    y <- cbind(
      c(y1, rep(NA, n)),  # Gaussian values in first half, NA in second
      c(rep(NA, n), y2)   # NA in first half, count values in second
    )
    
    # Data frame for covariates or index variables
    data <- data.frame(
      idx = c(Metadata_sub$idx1,Metadata_sub$idx1+max(Metadata_sub$idx1)),
      idx_actv = c(Metadata_sub$idx1,Metadata_sub$idx1),
      tcont = c(Metadata_sub$Group,rep(NA, n)),#,
      #x = factor(Metadata_sub[,Sample_id]),
      tactv = c(rep(NA, n),Metadata_sub$Group)
    )
    data$tcont = as.factor(data$tcont)
    data$tactv = as.factor(data$tactv)
    
    n_total <- nrow(Q_sparse)  # total length of the latent field
    
    # Create constraint matrix (each row is a constraint)
    
    A_constr <- Matrix(0, nrow = 2, ncol =  2*n_total, sparse = TRUE)
    A_constr[1,1:n_total] =1
    A_constr[2,(1:n_total)+n_total] =1
    
    # Right-hand side of constraints (all zero)
    e_constr <- c(0,0)
    
    #Q_inv <- MASS::ginv(as.matrix(Q_sparse))
    #avg_var <- mean(diag(Q_inv))
    #Q_sparse <- Q_sparse / avg_var
    
    
    aux_result <- try({
      
      res2 <- suppressWarnings(
        inla(
          y ~ 1+tcont +tactv+ f(idx, model = "generic0", Cmatrix = bdiag(Q_sparse,Q_sparse),
                         extraconstr = list(A = A_constr, e = e_constr)),
          data = data,
          family = c(countmodel , "binomial"),
          control.predictor = list(compute = TRUE,link = 1,link = 1),verbose = F,
          silent = TRUE,
          control.inla = list(control.vb = list(emergency = 500)),
          control.compute = list(dic = TRUE, waic = TRUE)
        )
      )
      
    }, silent = TRUE) 
    
    }else{
      
      
      y1 =  all_data[j,]
      y2 = ifelse(all_data[j,]==0,0,1)
      n = length(y1)
      # Response list: each element of the list is a vector of length n * 2 (for 2 likelihoods)
      if(countmodel=="binomial"){
        y <- y2
      }else{
        y <- y1+0.0001
      }
   
      # Data frame for covariates or index variables
      data <- data.frame(
        idx =  Metadata_sub$idx1,
        tcont = Metadata_sub$Group
      )
      data$tcont = as.factor(data$tcont)
      
      n_total <- nrow(Q_sparse)  # total length of the latent field
      
      # Create constraint matrix (each row is a constraint)
      
      A_constr <- Matrix(0, nrow = 1, ncol =  1*n_total, sparse = TRUE)
      A_constr[1,1:n_total] =1
      
      # Right-hand side of constraints (all zero)
      e_constr <- c(0)
    
      
      aux_result <- try({
        
        res2 <- suppressWarnings(
          inla(
            y ~ 1+tcont + f(idx, model = "generic0", Cmatrix = bdiag(Q_sparse),
                                  extraconstr = list(A = A_constr, e = e_constr)),
            data = data,
            family = c(countmodel),
            control.predictor = list(compute = TRUE,link = 1,link = 1),verbose = F,
            silent = TRUE,
            control.inla = list(control.vb = list(emergency = 500)),
            control.compute = list(dic = TRUE, waic = TRUE)
          )
        )
        
      }, silent = TRUE) 
      
      
 
    }
    
    if (class(aux_result)=="try-error") {
      message(paste("Error at iteration","gene",j))
      
      sig_cont[j] =NA
      sig_tactv[j] = NA
      treeEffect_cont[,j] = NA
      treeEffect_actv[,j] = NA
      SNR[j] =NA
      
      if(CollectPostD){
        
        PosteriorDist[[j]] =  NA
        
      }
      
      next
    }
    
    
    res2 = aux_result
    
    if(CollectPostD){
      
      PosteriorDist[[j]] =  res2
      
    }
    
    
    sig_cont[j] = !(res2$summary.fixed["tcont",3]<= 0 & res2$summary.fixed["tcont",5]>=0)
    sig_tactv[j] = !(res2$summary.fixed["tactv",3]<= 0 & res2$summary.fixed["tactv",5]>=0)
    
    beta_mean = res2$summary.fixed["tcont","mean"]
    beta_sd = res2$summary.fixed["tcont","sd"]
    z_value <- beta_mean / beta_sd
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    pvalue[1,j] = p_value
    
    beta_mean = res2$summary.fixed["tactv","mean"]
    beta_sd = res2$summary.fixed["tactv","sd"]
    z_value <- beta_mean / beta_sd
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    pvalue[2,j] = p_value
    
    treeEffect_cont[,j] = res2$summary.random$idx$mean[1:vcount(mst)]
    treeEffect_actv[,j] = res2$summary.random$idx$mean[1:vcount(mst)+max(1:vcount(mst))]
    SNR[1,j] = 1/res2$summary.hyperpar["Precision for idx","mean"]
    SNR[2,j] = 1/res2$summary.hyperpar["Precision for idx_actv","mean"]
  }
  adj_pvalue = pvalue
  adj_pvalue[1,] = p.adjust(pvalue[1,], method="BY")%>%round(digits = 3)
  adj_pvalue[2,] = p.adjust(pvalue[2,], method="BY")%>%round(digits = 3)
  
  return(list(sig_cont=sig_cont,
              sig_tactv = sig_tactv,
              treeEffect_cont=treeEffect_cont,
             treeEffect_actv = treeEffect_actv,
              SNR=SNR,
              mst=mst,
              Centers=Centers,
              Res = Res,
              pvalue =pvalue,
              adj_pvalue_BY = adj_pvalue,
              Metadata = Metadata_sub,
              PosteriorDist = PosteriorDist
  ))
}

Hcluster <- function(DataTocluster,thresholdGini=0.2,k,ClusterName){
  # Goal: Function to cluster a data using hierarchical clustering
  # INPUT: 
  # DataTocluster - > data matrix to cluster
  # thresholdGini - > Gini hyperparameter
  # OUTPUT: Cluste labled class
  
  cluster  <- genie::hclust2(objects=as.matrix(DataTocluster), thresholdGini=thresholdGini)
  clust    <- cutree(cluster, k = k)
  
  clust = as.matrix(clust) 
  colnames(clust) = ClusterName
  ClusteredData = cbind(clust,DataTocluster)
  
  return(ClusteredData)
}

AddColumnToFCS = function(InputDirOfData, ColumToadd,ColName,OutputDirOfData){
  # Goal: Function to add a column to an FCS file
  
  h <- spade:::SPADE.read.FCS(InputDirOfData)
  
  params <- flowCore::parameters(h)
  params1 <- params
  desc <- description(h)
  pd <- pData(params1)
  
  in_data <-  exprs(h)
  channel_number <- ncol(in_data) + 1
  channel_id <- paste("$P", channel_number, sep = "")
  channel_name <- ColName
  channel_range <- k + 1
  plist <- matrix(c(channel_name, channel_name, channel_range,
                    0, channel_range - 1))
  rownames(plist) <- c("name", "desc", "range", "minRange",
                       "maxRange")
  colnames(plist) <- c(channel_id)
  pd <- rbind(pd, t(plist))
  pData(params1) <- pd
  out_data <- cbind(in_data, ColumToadd)
  colnames(out_data) <- c(colnames(in_data),ColName)
  out_frame <- flowFrame(out_data, params1, description = desc)
  keyval <- list()
  keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
  keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
  keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
  keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
  keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
  keyword(out_frame) <- keyval
  write.FCS(out_frame, OutputDirOfData)
}


IntegrateTreeWithPriorTree <- function(ExprsData,ClusterCol,PriorColumn,  weighted.prob = 0.5,scale=0.1){
  
  # Determine the clusters with maximum similar cells between the data cluster
  # and manual clusters
  in_data   = data.frame(ExprsData)
  cluster   = ClusterCol
  spatialid = PriorColumn
  grp = c(cluster,spatialid)
  ######################
  in_data_grp = in_data %>% group_by_at(grp)
  in_data_grp_crossTab = in_data_grp %>% summarise(n = n()) %>% 
    ungroup() %>% group_by_at(cluster)%>%
    filter(n==max(n))
  
  n2=runif(nrow(in_data_grp_crossTab),0,1)
  in_data_grp_crossTab = in_data_grp_crossTab %>% ungroup%>%
    mutate(n2=n2)%>% group_by_at(cluster)%>%
    filter(n2==max(n2))
  
  
  # Compute the weighted average between centers of data and prior centers
  
  AdjustedCenters = matrix(0,nrow=nrow(in_data_grp_crossTab),ncol=ncol(in_data)-1)
  
  colnames(AdjustedCenters) = colnames(in_data[,!(colnames(in_data) %in% cluster)])
  
  auxA = in_data_grp_crossTab[,cluster] %>% as.matrix() %>% as.vector()
  auxB = in_data_grp_crossTab[,spatialid] %>% as.matrix() %>% as.vector()
  
  for(i in 1:nrow(in_data_grp_crossTab)){
    Aux_data_data = in_data %>% filter_at(vars(cluster),all_vars(.==auxA[i])) %>% select_at(vars(-cluster) )%>% mutate_at(vars(spatialid),function(x)x=0.001) %>% colMeans()
    Aux_data_prior = in_data %>% filter_at(vars(spatialid),all_vars(.==auxB[i])) %>% select_at(vars(-cluster) )%>% mutate_at(vars(spatialid),function(x)x=x*scale) %>% colMeans()
    
    aAux_data_prior = Aux_data_prior
    Aux_data_prior  = rep(0,length(aAux_data_prior))
    names(Aux_data_prior) = names(aAux_data_prior)
    Aux_data_prior[spatialid] =aAux_data_prior[spatialid]
    AdjustedCenters[i,] = weighted.prob * Aux_data_prior  +(1-weighted.prob)*Aux_data_data
    
  }
  return(list(UpdatedCenters=AdjustedCenters,ClusterAndPriorID = in_data_grp_crossTab))
}

ClusterToTree <- function(Centers,weighted = TRUE){
  # Derive a tree from a center
  # INPUT: Centers
  # OUTPUT: mst (igraph) object.
  adjacency <- dist(Centers, method = "manhattan")
  full_graph <- graph_from_adjacency_matrix(as.matrix(adjacency), mode = "undirected", 
                                weighted = weighted)
  mst <- minimum.spanning.tree(full_graph)
  
  return(mst)
}


ReconectDisconectedNetwk = function(mstdel){
  m2 = as_adjacency_matrix(mstdel)
  
  while(sum(rowSums(m2)==0)>0){
    
    Disc.Node =  which(rowSums(m2)==0)
    layoutCoord = layout_nicely(mstdel)
    dist = proxy::dist(layoutCoord[Disc.Node[1],,drop=F] ,layoutCoord)
    aux = which(as.vector(dist)==sort(as.vector(dist),decreasing =F)[2])
    mstdel= add_edges(mstdel, edges = c(Disc.Node[1], aux))
    m2 = as_adjacency_matrix(mstdel)
  }
  return(mstdel)
}

plotTree <- function(mst,PlotVar,
                     #community=FALSE,
                     main="Weighted",Lab=T,noLegend=F,
                     vertex.size=3,vertex.label.cex=1,
                     sizelegend="none",
                     limits = NULL,
                     cols=NULL,
                     legend.size = 0.2,
                     legend.text.size=10,
                     legend.size.alpha=1,
                     seed=11234){
  # Plot function
  # INPUT: mst -> graph object(igraph),community-> whether to plot community (logical),
  # main-> title,... Other parameters are igraph parameters. 
  # OUTPU: plot
  # Note, if PlotVar is a factor, plot is made from igraph base plot, if countitnous, ggplot is used.
  #if(community){
  #  a = cluster_fast_greedy(mst)
  #  V(mst)$id_com = membership(a)
  #  set.seed = seed
  #  plot(a,mst,vertex.size=8,vertex.label.cex=1,main=main)
  #}else{
  
  if(is.factor(PlotVar)){
    # a=colorRampPalette(c("blue","yellow","red"))
    #PlotVar = PlotVar %>%as.numeric()
    # aux_PlotVar =unique(PlotVar) %>% sort()
    # pal =  a(length(aux_PlotVar))[PlotVar]
    # set.seed = seed
    #  plot(mst,vertex.size=3,vertex.label.cex=1,vertex.label=NA,
    #       vertex.color =pal,main=main)
    #  legend("bottomright",bty = "n",
    #         legend=(aux_PlotVar),
    #        fill=(a(length(aux_PlotVar))), cex=.8,border=NA, horiz = TRUE,text.width = 0.001)
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link() + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_manual(values = cols)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(colour = guide_legend(override.aes = list(size=legend.size)),
             size = guide_legend(override.aes = list(alpha = 0.5))
      )
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }else{
    
    set.seed = seed
    pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link() + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_distiller(palette = "Spectral",limits=limits)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(#colour = guide_legend(override.aes = list(size=legend.size)),
             
             size = guide_legend(override.aes = list(alpha = alpha)))
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }
}


plotTree <- function(mst,PlotVar,
                     #community=FALSE,
                     main="Weighted",Lab=T,noLegend=F,
                     vertex.size=3,vertex.label.cex=1,
                     sizelegend="none",
                     edge_color="black",
                     edge_alpha=1,
                     limits = NULL,
                     cols=NULL,
                     legend.size = 0.2,
                     legend.text.size=10,
                     legend.size.alpha=1,
                     seed=11234){
  # Plot function
  # INPUT: mst -> graph object(igraph),community-> whether to plot community (logical),
  # main-> title,... Other parameters are igraph parameters. 
  # OUTPU: plot
  # Note, if PlotVar is a factor, plot is made from igraph base plot, if countitnous, ggplot is used.
  #if(community){
  #  a = cluster_fast_greedy(mst)
  #  V(mst)$id_com = membership(a)
  #  set.seed = seed
  #  plot(a,mst,vertex.size=8,vertex.label.cex=1,main=main)
  #}else{
  
  if(is.factor(PlotVar)){
    # a=colorRampPalette(c("blue","yellow","red"))
    #PlotVar = PlotVar %>%as.numeric()
    # aux_PlotVar =unique(PlotVar) %>% sort()
    # pal =  a(length(aux_PlotVar))[PlotVar]
    # set.seed = seed
    #  plot(mst,vertex.size=3,vertex.label.cex=1,vertex.label=NA,
    #       vertex.color =pal,main=main)
    #  legend("bottomright",bty = "n",
    #         legend=(aux_PlotVar),
    #        fill=(a(length(aux_PlotVar))), cex=.8,border=NA, horiz = TRUE,text.width = 0.001)
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link(color=edge_color,alpha=edge_alpha) + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_manual(values = cols)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(colour = guide_legend(override.aes = list(size=legend.size)),
             size = guide_legend(override.aes = list(alpha = 0.5))
      )
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }else{
    
    set.seed = seed
    pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))
    p= ggraph(mst, layout = "stress") + 
      geom_edge_link(color=edge_color,alpha=edge_alpha) + 
      geom_node_point(aes(color=PlotVar,size=vertex.size))+ guides(size=sizelegend)+
      scale_color_distiller(palette = "Spectral",limits=limits)+labs(color="",x="",y="",title=main,size="")+
      #  geom_node_text(label =V(mst),
      #                 colour = 'black', size=3,
      #                 show.legend = FALSE, family = "serif")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))+
      guides(#colour = guide_legend(override.aes = list(size=legend.size)),
        
        size = guide_legend(override.aes = list(alpha = alpha)))
    
    #   par(old.par)
    if(Lab){
      p = p +  geom_node_text(label =V(mst),
                              colour = 'black', size=3,
                              show.legend = FALSE, family = "serif")
    }
    if(noLegend){
      p=p+ guides(color="none")
    }
    if(sizelegend=="none"){
      p=p+ guides(size="none")
    }
    p
  }
}


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

c35 <- c(
  "dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "black",
  "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
  "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
  "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
  # Added 10 more
  "mediumpurple", "turquoise4", "darkseagreen", "tomato", "darkmagenta",
  "lightcoral", "midnightblue", "olivedrab", "sienna", "cyan4", "#E31A1C"
)
plotScatter = function(Xcord ,
                       Ycord,
                       Gene,
                       size=1,
                       main,
                       pal="Accent",
                       legend.size = 0.2,
                       legend.text.size=10,
                       noLegend=FALSE,
                       limits=NULL,
                       ManualColor =FALSE,
                       cols = NULL){
  
  if(!is.factor(Gene)){
    p = data.frame(Xcord,Ycord,Gene) %>% ggplot()+
      geom_point(aes(x=Xcord,y=Ycord,color=Gene),size=size)+
      scale_color_distiller(palette = "Spectral",limits=limits)+
      labs(color="",x="",y="",title=main)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            legend.text = element_text(size=legend.text.size))#+
    # +
    # guides(colour = guide_legend(override.aes = list(size=legend.size)))
  }else{
    
    p= data.frame(Xcord,Ycord,Gene) %>% ggplot()+
      geom_point(aes(x=Xcord,y=Ycord,color=Gene),size=size)+
      scale_color_brewer(palette=pal)+
      labs(color="",x="",y="",title=main)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            rect = element_blank(),
            legend.key.width= unit(0.2, 'cm'),
            #legend.key.size= unit(legend.size, 'cm'),
            legend.text = element_text(size=legend.text.size)
      )+
      guides(colour = guide_legend(override.aes = list(size=legend.size)))
  }
  if(noLegend){
    p=p+ guides(color="none")
  }
  if(ManualColor){
    p= p+scale_color_manual(values =cols )
  }
  p
}

#cv = function(v) {return((sd(v)*100)/mean(v))} 

ComputeCenter <- function(ExprsData,Genes,ClusterCol="cluster",f = mean){
  # Functions to compute centers
  # INPUT: ExprsData -> Expression data (cellxgenes), Genes-> Genes to which computation should be made
  #       ClusterCol -> Column with cluster labels, f-> Function to summarise for each cluster
  # OUTPU: Data frame with each Genes summarised with function f for each unique ClusterCol.
  
  ExprsData = as.data.frame(ExprsData)
  if(is.null(Genes)) Genes = colnames(ExprsData)[!(colnames(ExprsData)%in%ClusterCol)]
  G = c(ClusterCol,Genes)
  ExprsData =  ExprsData %>% select_at(G)
  aux =  ExprsData %>% group_by_at(ClusterCol) %>% summarise(
    across(where(is.numeric), list(f))
  )%>% ungroup
  colnames(aux) = ExprsData %>% ungroup %>% select_at(G) %>% colnames()
  
  return(aux)
}

IntegrateTwoTrees = function(ExprsDataRef,ExprsData,ClusterColRef,ClusterCol){
  
  # Funtion to update a Tree based on a new expression data
  # INPUT: ExprsDataRef-> Reference expression data(cellxgenes), ExprsData-> new expr data, 
  #        ClusterColRef -> cluster column for  ExprsDataRef, ClusterCol-> cluster column for ExprsData.
  # OUTPUT: Updated centers, outlier clusters id, and outlier id on updated clusters
  #Note: The two impute expressions must have the same genes;
  
  in_data_A = as.data.frame(ExprsDataRef)
  in_data_B = as.data.frame(ExprsData)
  oo = sort(unique(in_data_A[,ClusterColRef]))
  in_data_grp_crossTab = matrix(0, nrow = length(oo),ncol = 2)
  
  for (i in 1:length(oo)) {
    
    A =  in_data_A %>% dplyr::filter_at(vars(ClusterColRef),all_vars(.==oo[i])) %>% select_at(all_of(vars(-ClusterColRef))) %>% colMeans()
    #A$id = 1
    o = in_data_B[,ClusterCol] %>% unique
    pv=vector("numeric",length = length(o))
    
    sim = pv = NULL
    for (x in 1:length(o)){ 
      
      xy = o[x]
      B =  in_data_B %>% dplyr::filter_at(vars(ClusterCol),all_vars(.==xy)) %>% select_at(vars(-ClusterCol)) %>% colMeans()
      
      C   =  rbind(A,B)
      
      sim = c(sim, as.vector(proxy::simil(C, method = "manhattan")))
      
      
      
    }
    
    max.table2 = o[which(sim==max(sim,na.rm = T))]
    in_data_grp_crossTab[i,] =c(oo[i],max.table2)
    
  }
  
  colnames(in_data_grp_crossTab) = c("cluster_A","cluster_B")
  in_data_grp_crossTab = as.data.frame(in_data_grp_crossTab)
  
  # order according to A
  in_data_grp_crossTab = in_data_grp_crossTab[order(in_data_grp_crossTab$cluster_A),]
  
  # Determine clusters in B that are not captured
  outlier_cluster = unique(in_data_B[,ClusterCol])[!(unique(in_data_B[,ClusterCol]) %in% in_data_grp_crossTab$cluster_B)]
  outlier_center = in_data_B %>% dplyr::filter_at(vars(ClusterCol),all_vars(. %in% outlier_cluster)) %>% group_by_at(ClusterCol) %>%
    summarise(
      across(where(is.numeric), list(mean))
    ) %>% ungroup %>%  select_at(vars(-ClusterCol))
  colnames(outlier_center) = gsub("_1","",colnames(outlier_center))
  outlier_spatialid = max(in_data_A[,ClusterColRef])+(1:length(outlier_cluster))
  
  
  # Compute the weighted average between centers of data and prior ceenters
  weighted.prob = 0.5
  AdjustedCenters = matrix(0,nrow=nrow(in_data_grp_crossTab),ncol=ncol(in_data_A)-1)
  colnames(AdjustedCenters) = colnames(in_data_A %>% select_at(vars(-ClusterColRef)))
  for(i in 1:nrow(in_data_grp_crossTab)){
    Aux_data_data = in_data_A %>% dplyr::filter_at(vars(ClusterColRef),all_vars(.==in_data_grp_crossTab$cluster_A[i])) %>% dplyr::select_at(vars(-ClusterColRef)) %>% colMeans()
    Aux_data_prior = in_data_B %>% dplyr::filter_at(vars(ClusterCol),all_vars(.==in_data_grp_crossTab$cluster_B[i]))%>% dplyr::select_at(vars(-ClusterCol)) %>% colMeans()
    
    AdjustedCenters[i,] = weighted.prob * Aux_data_prior  + (1-weighted.prob) * Aux_data_data
    
  }
  AdjustedCenters = rbind(AdjustedCenters,outlier_center)
  
  return(list(UpdatedCenters=AdjustedCenters,
              unmachedCluster = outlier_cluster,
              IdOfUnmachedClusterOnUpdatedCenters=outlier_spatialid ))
}


IntegrateTwoTrees.impute = function(ExprsDataRef,ExprsData,ClusterColRef,ClusterCol){
  # Function to update a Tree based on a new expression data if some Genes are different.
  # INPUT: ExprsDataRef-> Reference expression data(cellxgenes), ExprsData-> new expr data, 
  #        ClusterColRef -> cluster column for  ExprsDataRef, ClusterCol-> cluster column for ExprsData.
  # OUTPUT: Updated centers, outlier clusters id, outlier id on updated clusters, and imputed expressions
  LabRef  = colnames(ExprsDataRef)
  Lab     = colnames(ExprsData)
  InT     = intersect(Lab,LabRef)
  FLab    = LabRef[!(LabRef%in%InT)]
  FLabRef = Lab[!(Lab%in%InT)]
  for (i in seq_len(length(FLabRef))) ExprsDataRef[,FLabRef[i]] = NA
  for (i in seq_len(length(FLab))) ExprsData[,FLab[i]] = NA
  
  in_data_A.mis  = as.data.frame(ExprsDataRef)
  in_data_B.mis  = as.data.frame(ExprsData)
  
  in_data_A.mis$id = 1
  in_data_B.mis$id = 2
  in_data.mis = rbind(in_data_A.mis,in_data_B.mis)
  
  
  # Remove cols that are all Missing (that is, the variable is missing in both data)
  resB = colSums(in_data.mis,na.rm = TRUE)
  in_data.mis = in_data.mis[,names(resB)[resB !=0]]
  
  # Set Missing col with missing as respons variable 
  
  resA = colSums(in_data.mis)
  aux.namesOfMis = names(resA[is.na(resA)])
  FiTTed = list()
  for(i in seq_len(length(aux.namesOfMis))){
    Y = in_data.mis %>% dplyr::select_at(vars(aux.namesOfMis[1]))
    Y = Y+1e-5 # Computational reasons
    X = in_data.mis %>% select (-all_of(aux.namesOfMis))
    
    DATA_imputation = cbind(Y,X)
    
    formula = paste(colnames(X%>%select(-id)%>% select_at((vars(-ClusterColRef,-ClusterColRef)))),
                    collapse = "+")
    formula = paste0(colnames(Y),"~",formula)
    formula = paste0(formula,"+f(id,model='iid')")
    formula = paste0(formula,"+f(",ClusterColRef,",model='iid')")
    if(ClusterColRef!=ClusterCol) formula = paste0(formula,"+f(",ClusterCol,",model='iid')")
    formula = as.formula(formula)
    
    Result = inla(formula = formula,data = DATA_imputation,family = "gamma",control.predictor=list(compute=TRUE),
    )
    aux = Result$summary.fitted.values$mean[is.na(Y)] %>% exp
    Y[is.na(Y)] = aux
    FiTTed[colnames(Y)][[1]] = aux
    
    ####################################### Using semi-structured modeling
    #formula = paste(colnames(X),
    #                collapse = "+")
    #formula = paste0("nn(~",formula,",size=",5,",decay=",0.01,")")
    #formula = paste0(colnames(Y),"~",formula)
    #formula = as.formula(formula)
    #Result = gamlss(formula, data=DATA_imputation,family ="GA")
    #Y[is.na(Y)] = fitted(Result)[is.na(Y)]
    ########################################
    in_data.mis[,colnames(Y)] = Y
    # Remove already imputed variable
    
    aux.namesOfMis = aux.namesOfMis[aux.namesOfMis!=colnames(Y)]
    
  }
  ##### End of imputation #########
  
  in_data_A = in_data.mis %>% filter(id==1) %>% select(-id)
  in_data_B = in_data.mis %>% filter(id==2) %>% select(-id)
  
  aux = IntegrateTwoTrees(ExprsDataRef =in_data_A ,ExprsData = in_data_B,
                          ClusterColRef ="clusterr" ,
                          ClusterCol = "clusterr")
  
  return(list(aux,Imputed = FiTTed))
  
}



GetTreeVariableGenes <- function(mst,
                                 ExprsData,
                                 ClusterCol, 
                                 useWeight=TRUE, 
                                 Robust = TRUE,
                                 Model="GA",
                                 nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - WeightGraphAdj
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the orriginal data to get graph info into the data
  id = data.frame(cluster = Idclust, 
                  Graphid = seq_len(length(Idclust)))
  id[,ClusterCol] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id,by = ClusterCol)
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% select(-Graphid) %>% colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted = in_data_AA[,names_gene]
  if(is.na(sum(in_data_AA))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = matrix(NA,ncol = length(names_gene))
  colnames(SNR) = names_gene
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene))
  colnames(BIC) = names_gene
  # Configure to hold tree-effect
  treeEffect = matrix(0,nrow =nrow(WeightGraphAdj) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% 
      dplyr::mutate(gene= gene+1e-02)
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    
    
    if (i==length(names_gene)) cat(' Finalizing...',"\n")
    if(!Robust) {
      Model=Model
    }else{
      m = fitDist(model_dataFrame$gene,type = "realplus")
      Model = m$family[1]
    }
    
    
    m1 <- gamlss(gene ~ gmrf(Graphid,
                             precision=as.matrix(PrecWeightGraphAdj)),
                 data=model_dataFrame,
                 family = Model)
    
    
    # Get Singnal to Noise ratio
    sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
    sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
    
    SNR[,names_gene[i]] = sSpat/sError
    BIC[,names_gene[i]] =  m1$sbc #-gen.likelihood(m1)()
    # if(i%%500 ==0) write.csv(SNR[,1:i] ,file=paste0(getwd(),"/SNR.txt"))
    
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = fitted(m1)
    # Get tree effects
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         fitted= in_data_AA.fitted[, names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # Get tree effects
    treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,BIC=BIC, Fitted =in_data_AA.fitted, TreeEffect = treeEffect ))
}


GetTreeVariableGenesLocScale <- function(mst,
                                         ExprsData,
                                         ClusterCol, 
                                         useWeight=TRUE, 
                                         Robust = TRUE,
                                         Model="GA",
                                         nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - WeightGraphAdj
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the orriginal data to get graph info into the data
  id = data.frame(cluster = Idclust, 
                  Graphid = seq_len(length(Idclust)))
  id[,ClusterCol] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id,by = ClusterCol)
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>% select(-Graphid) %>% colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted = in_data_AA[,names_gene]
  if(is.na(sum(in_data_AA))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = matrix(NA,ncol = length(names_gene))
  SNRScale = matrix(NA,ncol = length(names_gene))
  colnames(SNR) = names_gene
  colnames(SNRScale) = names_gene
  # Configure to hold tree-effect
  treeEffect = matrix(0,nrow =nrow(WeightGraphAdj) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% 
      dplyr::mutate(gene= gene+1e-02)
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    
    
    if (i==length(names_gene)) cat(' Finalizing...',"\n")
    if(!Robust) {
      Model=Model
    }else{
      m = fitDist(model_dataFrame$gene,type = "realplus")
      Model = m$family[1]
    }
    
    if("sigma" %in% m$parameters){
      m1 <- gamlss(gene ~ gmrf(Graphid,
                               precision=as.matrix(PrecWeightGraphAdj)),
                   ~   gmrf(Graphid,
                            precision=as.matrix(PrecWeightGraphAdj)),
                   data=model_dataFrame,
                   family = Model)
      
      # Get location Signal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[,names_gene[i]] = sSpat/sError
      
      #  Get Scale Signal to Noise ratio
      sError = coef(m1$sigma.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$sigma.coefSmo[[1]])[[2]] %>% exp
      
      SNRScale[,names_gene[i]] = sSpat/sError
      
    }else{
      
      m1 <- gamlss(gene ~ gmrf(Graphid,
                               precision=as.matrix(PrecWeightGraphAdj)),
                   data=model_dataFrame,
                   family = Model)
      
      # Get location Signal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[,names_gene[i]] = sSpat/sError
      
    }
    
    
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = fitted(m1)
    # Get tree effects
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         SNRScale=SNRScale[,names_gene[i]],
         fitted= in_data_AA.fitted[, names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    SNRScale[,names_gene[i]] = ko_param[[i]]$SNRScale
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # Get tree effects
    treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,SNRScale = SNRScale, Fitted =in_data_AA.fitted, TreeEffect = treeEffect ))
}


GetTreeVariableGenesDynamics <- function(mst,
                                         ExprsData,
                                         ClusterCol,
                                         TemporalCol,
                                         useWeight=TRUE, 
                                         Robust = TRUE,
                                         Model="GA",
                                         rho_tree = 0.9,
                                         rho_temp = 0.7,
                                         nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  
  idtime = sort(unique(ExprsData[,TemporalCol]))
  CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id      = id[-missId,]
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id, by = "TreeTemp")
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid) %>% 
    dplyr::select(-TreeTemp) %>% 
    colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted =in_data_AA.fitted2= in_data_AA[,names_gene,drop=F]
  if(is.na(sum(in_data_AA.fitted))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = SNR2 = matrix(NA,ncol = length(names_gene),nrow=2)
  colnames(SNR) = names_gene
  colnames(SNR2) = names_gene
  
  #Log-like 
  BIC = BIC2 = matrix(NA,ncol = length(names_gene))
  colnames(BIC) = names_gene
  colnames(BIC2) = names_gene
  # Configure to hold tree-effect
  treeEffect =treeEffect2 = matrix(NaN,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  colnames(treeEffect2) = names_gene
  rownames(treeEffect) = id$Graphid
  rownames(treeEffect2) = id$Graphid
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss",
                                    "Matrix"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    ## Implement Zero inflation here
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% dplyr::mutate(geneBin= as.numeric(gene>0) )
    
    model_dataFrameNonBinary = model_dataFrame %>% filter(geneBin==1)
    missId = which(!sort(unique(model_dataFrame$Graphid)) %in% unique(model_dataFrameNonBinary$Graphid))
    # id required to capture thevtree effect below
    auxid=as.numeric(as.character(sort(unique(model_dataFrameNonBinary$Graphid))))
    
    if(length(missId)>0) {
      TreeTemporaPrecision_count = TreeTemporaPrecision
      TreeTemporaPrecision_count = TreeTemporaPrecision_count[-missId,-missId]
      TreeTemporaPrecision_count = forceSymmetric(TreeTemporaPrecision_count)
    }
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    model_dataFrameNonBinary$Graphid=as.factor(model_dataFrameNonBinary$Graphid)
    
    #Check if zero inlated
    TruZero = 1-mean(model_dataFrame$geneBin)
    
    if(TruZero>0){
      if (i==length(names_gene)) cat(' Finalizing...',"\n")
      if(!Robust) {
        Model=Model
      }else{
        # m = fitDist(model_dataFrameNonBinary$gene,type = "realplus")
        m = fitDist(model_dataFrame$gene[model_dataFrame$gene>0],type = "realplus")
        Model = m$family[1]
      }
      
      
      #m1 <- gamlss(gene ~ gmrf(Graphid,
      #                         precision=as.matrix(TreeTemporaPrecision_count)),
      #             data=model_dataFrameNonBinary,
      #             family = Model)
      
      model_dataFrame$weight = 0
      model_dataFrame$weight[model_dataFrame$geneBin==0] = (1/sum(model_dataFrame$geneBin==0))*1e+3
      model_dataFrame$weight[model_dataFrame$geneBin==1] = 1/sum(model_dataFrame$geneBin==1)*1e+3
      model_dataFrame$gene = model_dataFrame$gene+1e-2
      model_dataFrameNonBinary$weight =  1/sum(model_dataFrame$geneBin==1)*1e+3
      
      #  m2 <- gamlss(geneBin ~ gmrf(Graphid,
      #                           precision=as.matrix(TreeTemporaPrecision)),
      #               data=model_dataFrame,weights = weight,
      #               family = "BI")
      
      
      
      m1 <-     gamlss(gene ~ gmrf(Graphid,
                                   precision=as.matrix(TreeTemporaPrecision_count)),
                       data=model_dataFrameNonBinary,#weights = weight,
                       family = Model)
      
      # Get Signal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[1,names_gene[i]] = sSpat
      SNR[2,names_gene[i]] = sError
      
      BIC[,names_gene[i]] =  m1$sbc 
      
      # sError = coef(m2$mu.coefSmo[[1]])[[1]] %>% exp
      # sSpat = coef(m2$mu.coefSmo[[1]])[[2]] %>% exp
      
      # SNR2[1,names_gene[i]] = sSpat/sError
      #  BIC2[,names_gene[i]] =  m1$sbc 
      #-gen.likelihood(m1)()
      # if(i%%500 ==0) write.csv(SNR[,1:i] ,file=paste0(getwd(),"/SNR.txt"))
      
      # Get smoothed/predicted gene expression
      #in_data_AA.fitted[model_dataFrame$geneBin==1, names_gene[i]] = fitted(m1)
      #  in_data_AA.fitted[, names_gene[i]] = fitted(m1)
      #  in_data_AA.fitted2[, names_gene[i]] = fitted(m2)
      # Get tree effects
      
      #treeEffect[as.character(auxid),names_gene[i]]  = getSmo(m1)$beta
      #  treeEffect[ ,names_gene[i]]  = getSmo(m1)$beta
      #  treeEffect2[,names_gene[i]]  = getSmo(m2)$beta
      
    }else{
      model_dataFrame$gene = model_dataFrame$gene+1e-3
      
      if(!Robust) {
        Model=Model
      }else{
        m = fitDist(model_dataFrame$gene,type = "realplus")
        Model = m$family[1]
      }
      
      
      m1 <- gamlss(gene ~ gmrf(Graphid,
                               precision=as.matrix(TreeTemporaPrecision)),
                   data=model_dataFrame,
                   family = Model)
      
      # Get Singnal to Noise ratio
      sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
      sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
      
      SNR[,names_gene[i]] = sSpat/sError
      BIC[,names_gene[i]] =  m1$sbc 
      
      
      
      # Get smoothed/predicted gene expression
      # in_data_AA.fitted[, names_gene[i]] = fitted(m1)
      # Get tree effects
      
      # treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    }
    
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],SNR2=SNR2[,names_gene[i]],
         BIC = BIC[,names_gene[i]],BIC2 = BIC2[,names_gene[i]]
         #  fitted= in_data_AA.fitted[, names_gene[i]],fitted2= in_data_AA.fitted2[, names_gene[i]],
         #  treeEffect = treeEffect[,names_gene[i]],treeEffect2 = treeEffect2[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    SNR2[,names_gene[i]] = ko_param[[i]]$SNR2
    BIC2[,names_gene[i]] = ko_param[[i]]$BIC2
    # Get smoothed/predicted gene expression
    #in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # in_data_AA.fitted2[, names_gene[i]] = ko_param[[i]]$fitted2
    # Get tree effects
    # treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
    # treeEffect2[,names_gene[i]]  = ko_param[[i]]$treeEffect2
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,SNR2 = SNR2,
              BIC=BIC, BIC2=BIC2))
  # Fitted =in_data_AA.fitted, Fitted2 =in_data_AA.fitted2,
  # TreeEffect = treeEffect,TreeEffect2 = treeEffect2 ))
}



GetTreeVariableGenesDynamicsMisHandler <- function(mst,
                                                   ExprsData,
                                                   ClusterCol,
                                                   TemporalCol,
                                                   useWeight=TRUE, 
                                                   Robust = TRUE,
                                                   Model="GA",
                                                   rho_tree = 0.9,
                                                   rho_temp = 0.7,
                                                   nCores=1){
  
  
  Idclust = sort(unique(ExprsData[,ClusterCol]))
  nclust = length(Idclust)
  nNodes = vcount(mst)
  
  if(nclust!=nNodes) stop("Nodes on data and graph are different")
  if(useWeight){  WeightGraphAdj = as_adjacency_matrix(mst,attr = "weight")
  }else{  WeightGraphAdj = as_adjacency_matrix(mst)}
  
  rho_tree = 0.9
  PrecWeightGraphAdj = Diagonal(nrow(WeightGraphAdj),x=rowSums(WeightGraphAdj!=0)) - rho_tree * WeightGraphAdj
  
  ## Temporal
  rho_temp = 0.7
  idtime = sort(unique(ExprsData[,TemporalCol]))
  CovTemporal = Q.AR1(length(idtime),1,rho_temp,vcov = T)
  
  TreeTemporaCov = kronecker(solve(PrecWeightGraphAdj),CovTemporal)
  TreeTemporaPrecision = solve(TreeTemporaCov) 
  
  # Tree & Temporal
  ExprsData[,"TreeTemp"] = (paste0(ExprsData[,ClusterCol],ExprsData[,TemporalCol]))
  
  idTreeTime = expand.grid(Idclust,idtime)
  idTreeTime = paste0(idTreeTime$Var1,idTreeTime$Var2)
  
  # Since the mst_graph_fromPriorKwldge was sorted according to in_data_grp_crossTab$cluster_A
  # We can link it to the original data to get graph info into the data
  id = data.frame(cluster = idTreeTime, 
                  Graphid = seq_len(length(idTreeTime)))
  
  missId  = id$Graphid[! id$cluster %in% unique(ExprsData[,"TreeTemp"])]
  id      = id[-missId,]
  id[,"TreeTemp"] = id[,"cluster"]
  id = id %>% select(-cluster)
  
  TreeTemporaPrecision = TreeTemporaPrecision[-missId,-missId]
  TreeTemporaPrecision  = forceSymmetric(TreeTemporaPrecision)
  # Update input data
  
  in_data_AA = ExprsData  %>% left_join(id, by = "TreeTemp")
  
  # Calculate gene signal to noise ratio
  
  names_gene = in_data_AA %>% dplyr::select(-all_of(ClusterCol)) %>%
    dplyr::select(-all_of(TemporalCol))  %>% 
    dplyr::select(-Graphid) %>% 
    dplyr::select(-TreeTemp) %>% 
    colnames()
  
  # Configure to hold predicted gene values
  in_data_AA.fitted = in_data_AA[,names_gene]
  if(is.na(sum(in_data_AA.fitted))) warning("Expression contains missing data")
  
  # Configured to hold Signal-to-noise ratio
  SNR = matrix(NA,ncol = length(names_gene))
  colnames(SNR) = names_gene
  
  #Log-like 
  BIC = matrix(NA,ncol = length(names_gene))
  colnames(BIC) = names_gene
  # Configure to hold tree-effect
  treeEffect = matrix(0,nrow =nrow(TreeTemporaPrecision) ,ncol = length(names_gene))
  colnames(treeEffect) = names_gene
  rownames(treeEffect) = id$Graphid
  # Compute density to derive nodes weights
  
  # weight = in_data_AA %>% group_by_at(ClusterCol) %>% dplyr::summarise(n=n()) %>%
  #          group_by_at(ClusterCol) %>%
  #         dplyr::reframe(dense = n/nrow(in_data_AA),weight = 1/dense) %>% select(-dense)
  # weight$weight = weight$weight/sum(weight$weight)
  #weight = 1/dense
  #weight = weight/sum(weight)
  #  in_data_AA = in_data_AA%>% left_join(weight,by=ClusterCol)
  # Parallel computing
  
  
  cl <- makeCluster(nCores, outfile="")
  registerDoParallel(cl)
  
  
  ko_param <- foreach(i=seq_len(length(names_gene)),.errorhandling = "pass",
                      .packages = c("doParallel",
                                    "Rfast2",
                                    "tidyverse",
                                    "gamlss.spatial",
                                    "gamlss"
                      )
  ) %dopar% {
    
    cat("Modeling gene", (i/length(names_gene))*100,"%","\n" )
    
    
    model_dataFrame = in_data_AA %>% 
      dplyr::select(names_gene[i], Graphid) %>% 
      dplyr::rename(gene =names_gene[i]) %>% 
      dplyr::mutate(gene= gene+1e-02)
    
    if (sum(is.na(model_dataFrame$gene)) > 0 | sum(is.nan(model_dataFrame$gene))>0 ){
      warning("NA in Gene",names_gene[i], ". Hence will be neglected") 
      SNR[,names_gene[i]] = NA
      in_data_AA.fitted[, names_gene[i]] = NA
      treeEffect[,names_gene[i]]  = NA
      next
    }
    
    model_dataFrame$Graphid = as.factor(model_dataFrame$Graphid)
    
    
    if (i==length(names_gene)) cat(' Finalizing...',"\n")
    if(!Robust) {
      Model=Model
    }else{
      m = fitDist(model_dataFrame$gene,type = "realplus")
      Model = m$family[1]
    }
    
    
    for(k in 1: length(names(m$fits))){
      
      Model =  names(m$fits)[k]
      m1 =  tryCatch(
        expr = {
          m1 <- gamlss(gene ~ gmrf(Graphid,
                                   precision=as.matrix(TreeTemporaPrecision)),
                       data=model_dataFrame,
                       family = Model)
        },
        error = function(e){ 
          m1 = NULL
        },
        warning = function(w){
          
          m1 = m1
          
        },
        finally = {
          #m1
        }
      )
      
      if(!is.null(m1)){
        break
      }
    }
    # Get Singnal to Noise ratio
    sError = coef(m1$mu.coefSmo[[1]])[[1]] %>% exp
    sSpat = coef(m1$mu.coefSmo[[1]])[[2]] %>% exp
    
    SNR[,names_gene[i]] = sSpat/sError
    BIC[,names_gene[i]] =  m1$sbc #-gen.likelihood(m1)()
    # if(i%%500 ==0) write.csv(SNR[,1:i] ,file=paste0(getwd(),"/SNR.txt"))
    
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = fitted(m1)
    # Get tree effects
    treeEffect[,names_gene[i]]  = getSmo(m1)$beta
    
    list(Gene= names_gene[i],
         SNR=SNR[,names_gene[i]],
         BIC = BIC[,names_gene[i]],
         fitted= in_data_AA.fitted[, names_gene[i]],
         treeEffect = treeEffect[,names_gene[i]]
    )
  }
  
  
  for (i in seq_len(length(names_gene))) {
    if(is.null(ko_param[[i]]$SNR)) next
    # Get smoothed/predicted gene expression
    SNR[,names_gene[i]] = ko_param[[i]]$SNR
    BIC[,names_gene[i]] = ko_param[[i]]$BIC
    # Get smoothed/predicted gene expression
    in_data_AA.fitted[, names_gene[i]] = ko_param[[i]]$fitted
    # Get tree effects
    treeEffect[,names_gene[i]]  = ko_param[[i]]$treeEffect
  }
  
  parallel::stopCluster(cl) 
  return(list(SNR = SNR,BIC=BIC, Fitted =in_data_AA.fitted, TreeEffect = treeEffect ))
}

# Dual Mesh

book.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

