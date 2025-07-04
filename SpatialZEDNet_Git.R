SpatialDENet <- function(all_data,
                         metadata,
                         Sample_id = "Sample_id",
                         CollectPostD = FALSE){
  library(INLA)
  library(sf)
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
    
    mesh = inla.mesh.2d(loc.domain = coord, offset = c(1, 1), 
                        max.edge = c(.7, .7), cutoff =1)
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
  treeEffect_actv = treeEffect_cont = matrix(NA, nrow =vcount(mst),ncol = n_genes)
  PosteriorDist = list()
  
  pb <- progress_bar$new(total = n_genes)
  
  for (j in 1:n_genes) {
    pb$tick()
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
      idx = c(Metadata_sub$idx1,Metadata_sub$idx1),
      tcont = c(Metadata_sub$Group,Metadata_sub$Group),
      x = factor(Metadata_sub[,Sample_id])
      #tactv = c(rep(NA, n),metadata$Group)
    )
    data$tcont = as.factor(data$tcont)
    #constrained_indices <- c(1:(2*dim(Q)[1]))
    n_total <- nrow(Q_sparse)  # total length of the latent field
    
    # Create constraint matrix (each row is a constraint)
    #A_constr <- Matrix(0, nrow = length(constrained_indices), ncol = n_total, sparse = TRUE)
    #for (i in seq_along(constrained_indices)) {
    #  A_constr[i, constrained_indices[i]] <- 1
    #}
    A_constr <- Matrix(0, nrow = 1, ncol =  n_total, sparse = TRUE)
    A_constr[1,1:n_total] =1
    #A_constr[2,(1:n_total)+n_total] =1
    # Right-hand side of constraints (all zero)
    e_constr <- c(0)
    
    aux_result <- try({
      
      res2 <- suppressWarnings(
        inla(
          y ~ 1+tcont + f(idx, model = "generic0", Cmatrix = bdiag(Q_sparse),
                          extraconstr = list(A = A_constr, e = e_constr)),
          data = data,
          family = c("lognormal", "binomial"),
          control.predictor = list(compute = TRUE,link = 1,link = 1),verbose = F,
          silent = TRUE,
          control.inla = list(control.vb = list(emergency = 500))
        )
      )
      
    }, silent = TRUE) 
    
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
    #sig_tactv[j] = !(res2$summary.fixed["tactv",3]<= 0 & res2$summary.fixed["tactv",5]>=0)
    
    beta_mean = res2$summary.fixed["tcont","mean"]
    beta_sd = res2$summary.fixed["tcont","sd"]
    
    z_value <- beta_mean / beta_sd
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    pvalue[j] = p_value
    
    treeEffect_cont[,j] = res2$summary.random$idx$mean[1:vcount(mst)]
    #treeEffect_actv[,j] = res2$summary.random$idx$mean[(1:vcount(mst))+vcount(mst)]
    SNR[j] = 1/res2$summary.hyperpar["Precision for idx","mean"]
    
  }
  
  
  return(list(sig_cont=sig_cont,
              # sig_tactv = sig_tactv,
              treeEffect_cont=treeEffect_cont,
              # treeEffect_actv = treeEffect_actv,
              SNR=SNR,
              mst=mst,
              Centers=Centers,
              Res = Res,
              mnean =pvalue,
              Metadata = Metadata_sub
  ))
}