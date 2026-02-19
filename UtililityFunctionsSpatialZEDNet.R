# Function Documentation
# Optimized functions
########################
# Packages Required
########################

SpatialDENet <- function(all_data,
                         metadata,
                         Sample_id = "Sample_id",
                         joint=TRUE,
                         countmodel = "nbinomial",
                         interface="INLA",
                         CollectPostD = FALSE,
                         offset = c(.2, .2),
                         max.edge = c(.6, .6),
                         cutoff =.1){
  
  
  
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
    aux = lapply(1:nrow(coordsPoints),function(x)st_intersects(coordsPoints[x,], dmesh2)[[1]][1])
    aux = aux %>% unlist()
    
    metadata_sub <- bind_cols(metadata_sub,meshid=paste0(aux,k))
    Metadata_sub = rbind(Metadata_sub,metadata_sub)
    
  }
  all_data = as.matrix(All_data)
  Metadata_sub$meshid_id = as.numeric(as.factor(Metadata_sub$meshid))
  coords = Metadata_sub[,c("x","y")]
  
  PC =  prcomp(t(all_data), center = TRUE, scale. = TRUE)$x
  
  Centers = bind_cols(col=Metadata_sub$meshid_id,PC[,1:5]) %>% group_by(col) %>% summarise_all(mean)%>%dplyr::select(-col)
  Res = data.frame(col=Metadata_sub$meshid_id) %>% group_by(col) %>%summarise(nn = n())
  
  mst = ClusterToTree(Centers)
  
  
  
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
  
  
  
  if(! interface=="IWLS"){
    
    ##### Model gene-by-gene
    
    pvalue = sig_cont = sig_tactv = SNR= rep(0, n_genes)
    names(pvalue) = names(sig_cont) = names(sig_tactv) = names(SNR) = rownames(all_data)
     SNR = matrix(NA,2,n_genes)
    pvalue =  matrix(NA,3,n_genes)
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
          idx = c(Metadata_sub$idx1,Metadata_sub$idx1),#+max(Metadata_sub$idx1)),
          idx_actv = c(Metadata_sub$idx1,Metadata_sub$idx1),
          tcont = c(Metadata_sub$Group,Metadata_sub$Group),#,
          #x = factor(Metadata_sub[,Sample_id]),
          tactv = c(rep(NA, n),Metadata_sub$Group)
        )
        data$tcont = as.factor(data$tcont)
        data$tactv = as.factor(data$tactv)
        
        n_total <- nrow(Q_sparse)  # total length of the latent field
        
        # Create constraint matrix (each row is a constraint)
        
        A_constr <- Matrix(0, nrow = 1, ncol =   n_total, sparse = TRUE)
        A_constr[1,1:n_total] =1
        # A_constr[2,(1:n_total)+n_total] =1
        
        # Right-hand side of constraints (all zero)
        e_constr <- c(0)
        
        #Q_inv <- MASS::ginv(as.matrix(Q_sparse))
        #avg_var <- mean(diag(Q_inv))
        #Q_sparse <- Q_sparse / avg_var
        
        
        aux_result <- try({
          
          res2 <- suppressWarnings(
            inla(
              y ~ 1+tcont +tactv+ f(idx, model = "generic0", Cmatrix = bdiag(Q_sparse),
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
        }else if(countmodel=="nbinomial"){
          y <- y1
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
              y ~ 1+tcont + f(idx, model = "generic1", Cmatrix = bdiag(Q_sparse),
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
      pvalue1 = p_value
      
      beta_mean = res2$summary.fixed["tactv","mean"]
      beta_sd = res2$summary.fixed["tactv","sd"]
      z_value <- beta_mean / beta_sd
      p_value <- 2 * (1 - pnorm(abs(z_value)))
      pvalue2 = p_value
      
      
      pval = c(pvalue1,pvalue2)
      
      pvalue[1,j] = pvalue1
      pvalue[2,j] = pvalue2
      pvalue[3,j] = cauchy_combine(pval[!is.na(pval)])
      
      treeEffect_cont[,j] = res2$summary.random$idx$mean[1:vcount(mst)]
      treeEffect_actv[,j] = res2$summary.random$idx$mean[1:vcount(mst)+max(1:vcount(mst))]
      SNR[1,j] = 1/res2$summary.hyperpar["Precision for idx","mean"]
      SNR[2,j] = 1/res2$summary.hyperpar["Precision for idx_actv","mean"]
    }
    adj_pvalue = pvalue
    adj_pvalue[1,] = p.adjust(pvalue[1,], method="BY")%>%round(digits = 3)
    adj_pvalue[2,] = p.adjust(pvalue[2,], method="BY")%>%round(digits = 3)
    adj_pvalue[3,] = p.adjust(pvalue[3,], method="BY")%>%round(digits = 3)
    
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
  }else if(interface=="IWLS" & countmodel == "binomial"){
    
    
    counts = all_data
    sampleid = metadata$Group
    spatialCov = Metadata_sub$idx1
    G <- nrow(counts); n <- ncol(counts)
    
    # computeSizeFactors <- function(countData) {
    # Step 1: compute geometric mean per gene
    #     geom_means <- apply(countData, 1, function(gene_counts) {
    # Avoid zeros by returning NA if all counts are zero
    #       if (all(gene_counts == 0)) return(NA)
    #        exp(mean(log(gene_counts[gene_counts > 0])))
    #      })
    
    # Step 2: compute ratios K_ij / g_i for each gene
    #     ratios <- countData / geom_means
    
    # Step 3: take median of ratios for each sample
    #     size_factors <- apply(ratios, 2, function(x) median(x[!is.na(x) & is.finite(x)], na.rm = TRUE))
    
    # Normalize size factors to have geometric mean = 1 (optional)
    #    size_factors <- size_factors / exp(mean(log(size_factors+0.0001)))
    
    #    return(size_factors)
    #  }
    
    # 1. size factors
    #s_j <- computeSizeFactors(counts)  # vector length n
    
    coldata = data.frame(condition = metadata$Group)
    # 2. design matrices
    cond <- as.factor(as.numeric(factor(coldata$condition-min(coldata$condition)))-1)
    X_full  <- model.matrix(~ cond)         # includes intercept + condition
    X_null  <- X_full[,1,drop=F]            # intercept-only model
    
    # 3. dispersion estimate per gene (simple method of moments; DESeq2 is preferable)
    # method-of-moments: alpha = max( (var - mean) / mean^2, small_eps )
    # estimate_alpha_mom <- function(k, s) {
    #   mu <- as.numeric(rowMeans(matrix(k / s, nrow = 1))) * s # simplified for 1 gene
    # better: compute normalized mean: g = mean(k/s)
    #    g <- mean(k / s)
    #     m <- g * mean(s) / mean(s) # reduce to g
    # simpler, use sample mean and var on counts:
    #     m0 <- mean(k)
    #     v0 <- var(k)
    #     alpha <- (v0 - m0) / (m0^2)
    #     if (is.na(alpha) || alpha < 1e-8) alpha <- 1e-8
    #     alpha
    #   }
    
    #   dds <- DESeqDataSetFromMatrix(countData = counts+1,
    #                                  colData = metadata[,"Group",drop=F],
    #                                  design = ~ Group)
    #    dds <- estimateSizeFactors(dds)
    ##    dds <- estimateDispersions(dds, fitType = "parametric")
    #    Alph = dispersions(dds)
    
    M = NULL
    for (i in 1:G) {
      K = counts[i,]
       #m =   glm.nb(K~.,data = data.frame(K=K,x=X_full[,"cond1"]))$coefficients["x"][[1]]
        m =   lm(K~.,data = data.frame(K=K,x=X_full[,"cond1"]))$coefficients["x"][[1]]
      M = c(M,m)
    }
    M = as.matrix(M)
    
    sigma_prior <- apply(M, 2, function(b) {
      quantile(abs(b), prob=0.95, na.rm = TRUE) / qnorm(0.975)
    })
    sigma2_prior <- sigma_prior
    
    mean_counts <- rowMeans(counts)
    beta_abs <- abs(M[,1])  # e.g. condition coefficient
    fit_trend <- loess(beta_abs ~ log1p(mean_counts), span = 0.5)
    sigma_prior_trend <- predict(fit_trend)
    sigma2_prior <- sigma_prior_trend^2
    
    # 4. function to compute negative binomial log-likelihood given fitted beta and alpha
    # nb_loglik <- function(K, mu, alpha) {
    #   # R's rnbinom uses 'size' = 1/alpha, with mean=mu
    #   size <- ifelse(alpha <= 0, 1e8, 1/alpha)
    #   sum(dnbinom(K, size=size, mu=mu, log=TRUE))
    #  }
    
    # BINOMIAL COMPONENT
    ZS = as.matrix(fastDummies::dummy_cols(spatialCov,remove_selected_columns = TRUE))
    
    # 5. Per-gene testing (parallel)
    ncores <- max(1, detectCores()-1)
    cl <- makeCluster(ncores)
    clusterExport(cl, c("counts","coords","X_full","X_null","ZS","Q_sparse",
                        "sigma2_prior","sampleid","fit_binorm_glm_spatial_R"))
    clusterEvalQ(cl, c(library(MASS),library(tidyverse),library(Matrix),library(sf),library(gstat))) # if used inside
    
    test_gene <- function(i) {
      
      sigma2 <- c(1e6, sigma2_prior[i])
      
      K = as.numeric((counts[i,]-1) >0)
      
     
      fit_null <-  fit_binorm_glm_spatial_R(K=K, X=X_null,sigma2=sigma2,
                                            coords =coords ,Q_sparse =Q_sparse ,ZS =ZS ,
                                            tau2 = c(.1,.1), phi = c(.1,.1),sampleid,spatial=T,
                                            intercept_unpenalized = TRUE,
                                            maxit = 100, tol = 1e-8, verbose = FALSE)
      
      

      mu_null <- as.numeric(1/(1+ exp(-(X_null %*% fit_null$beta+ZS%*%fit_null$u))))
      ll_null <- sum(dbinom(K , p = mu_null,size=1,log = T))
      
      # fit full
      
      fit_full <- fit_binorm_glm_spatial_R(K = K, X=X_full,sigma2=sigma2,
                                           coords,Q_sparse =Q_sparse ,ZS =ZS,
                                           tau2 = c(.01,.1), phi = c(1,1),sampleid,spatial=T,
                                           intercept_unpenalized = TRUE,
                                           maxit = 100, tol = 1e-8, verbose = FALSE)
      
      
      mu_full <- as.numeric(1/(1+ exp(-(X_full %*% fit_full$beta+ ZS %*%fit_full$u))))
      ll_full <- sum(dbinom(K , p = mu_full,size=1,log = T))
      
      
      D <- 2 * (ll_full - ll_null)
      pval <- pchisq(D, df = ncol(X_full) - ncol(X_null), lower.tail = FALSE)
      
      wald_pvalue = 2 * pnorm(-abs(fit_full$beta[2]/sqrt(fit_full$cov[2,2])))
      
      c(stat = D, p = pval,gene =rownames(counts[i,,drop=F]) )
    }
    
    res_list <- parLapply(cl, 1:G, test_gene)
    stopCluster(cl)
    
    res_mat <- do.call(rbind, res_list)
    gene = res_mat[, "gene"]
    res_mat <- lapply(1:2,function(i) res_mat[,i]%>%as.numeric) 
    res_mat <- do.call(cbind, res_mat)
    
    pvals <- res_mat[, 2]
    p.adj <- p.adjust(pvals, method = "BH")
    
    res_df <- data.frame(gene = gene, stat = res_mat[,1], pvalue = pvals, padj = p.adj)
    
    return(res_df)
    
  }else if(interface=="IWLS" & countmodel == "nbinomial"){
    
    counts = all_data
    sampleid = metadata$Group
    spatialCov = Metadata_sub$idx1
    G <- nrow(counts); n <- ncol(counts)
    
    # 1. size factors
    s_j <- computeSizeFactors(counts)  # vector length n
    
    coldata = data.frame(condition =metadata$Group)
    # 2. design matrices
    cond <- as.factor(as.numeric(as.factor(coldata$condition))-1)
    sampleid = as.factor(as.numeric(as.factor(coldata$condition))-1)
    X_full  <- model.matrix(~ cond)         # includes intercept + condition
    X_null  <- X_full[,1,drop=F]            # intercept-only model
    
    # 3. dispersion estimate per gene (simple method of moments; DESeq2 is preferable)
    # method-of-moments: alpha = max( (var - mean) / mean^2, small_eps )
    estimate_alpha_mom <- function(k, s) {
      mu <- as.numeric(rowMeans(matrix(k / s, nrow = 1))) * s # simplified for 1 gene
      # better: compute normalized mean: g = mean(k/s)
      g <- mean(k / s)
      m <- g * mean(s) / mean(s) # reduce to g
      # simpler, use sample mean and var on counts:
      m0 <- mean(k)
      v0 <- var(k)
      alpha <- (v0 - m0) / (m0^2)
      if (is.na(alpha) || alpha < 1e-8) alpha <- 1e-8
      alpha
    }
    
    
    ##### Estimate Alpha
    
    out_mom <- estimate_alpha_mom_vec(counts, size_factors = s_j)
    mu_bar <- out_mom$mu       # mean(normalized counts) per gene
    alpha_gene <- out_mom$alpha
    
    # 2. fit trend
    trend <- fit_dispersion_trend(mu = mu_bar, alpha = alpha_gene, use_loess = T)
    alpha_trend <- trend$fitted
    
    # 3. shrink
    n_samples <- ncol(counts)
    shrink_out <- shrink_dispersions(mu = mu_bar, alpha_gene = alpha_gene, alpha_trend = alpha_trend,
                                     n_samples = n_samples, n_params = ncol(X_full), sampling_var_method = "analytic")
    Alph <- shrink_out$alpha_final
    
    #####
    M = NULL
    for (i in 1:G) {
      K = counts[i,]
      #m =   glm.nb(K~.,data = data.frame(K=K, x=X_full[,"cond1"]))
      m =   lm(K~.,data = data.frame(K=K, x=X_full[,"cond1"]))
      m1 = m$coefficients["x"][[1]]
      M = c(M,m1)
    }
    M = as.matrix(M)
    
    sigma_prior <- apply(M, 2, function(b) {
      quantile(abs(b), prob=0.95, na.rm = TRUE) / qnorm(0.975)
    })
    sigma2_prior <- sigma_prior
    
    mean_counts1 <- rowMeans(counts[,sampleid==0])
    mean_counts2 <- rowMeans(counts[,sampleid==1])
    
    # 
    beta_abs <- abs(M[,1])  # e.g. condition coefficient
    fit_trend <- loess(beta_abs ~ log1p(mean_counts1)+log1p(mean_counts2), span = 0.5)
    
    #fit_trend <- loess(beta_abs ~ log1p(mean_counts1), span = 0.5)
    sigma_prior_trend <- predict(fit_trend)
    sigma2_prior <- sigma_prior_trend^2
    
    # 4. function to compute negative binomial log-likelihood given fitted beta and alpha
    nb_loglik <- function(K, mu, alpha) {
      # R's rnbinom uses 'size' = 1/alpha, with mean=mu
      size <- ifelse(alpha <= 0, 1e8, 1/alpha)
      sum(dnbinom(K, size=size, mu=mu, log=TRUE))
    }
    
    # NEGATIVE BINOMIAL
    
    ZS = as.matrix(fastDummies::dummy_cols(spatialCov,remove_selected_columns = TRUE))
    
    ncores <- max(1, detectCores()-1)
    cl <- makeCluster(ncores)
    clusterExport(cl, c("counts","coords","s_j","X_full","X_null","ZS","Q_sparse",
                        "fit_nb_glm_spatial_R","nb_loglik","estimate_alpha_mom","Alph","RobustTauPhi","sigma2_prior","sampleid","RobustTauPhi_Matern","cond"))
    clusterEvalQ(cl, c(library(MASS),library(tidyverse),library(Matrix),library(sf),library(gstat))) # if used inside
    
    test_gene <- function(i) {
      
      X_full  <- model.matrix(~ cond)         # includes intercept + condition
      X_null  <- X_full[,1,drop=F]            # intercept-only model
      
      K <- as.numeric(counts[i, ])
      #alpha = Alph[i]
      alpha <- estimate_alpha_mom(K, s_j)  # or use reliable estimator / DESeq2
      alpha = Alph[i]
      # choose prior variances: very weak prior on beta, spatial prior set by tau2 & phi
      sigma2 <- c(1e6, sigma2_prior[i])  # large => effectively unpenalized fixed effects
      # fit null
      
   
      
      fit_full <- fit_nb_glm_spatial_R(K = K, X = X_full, s = s_j, alpha = alpha,
                                       sigma2 = sigma2, coords = coords,Q_sparse=Q_sparse,ZS=ZS,
                                       tau2 = c(2,2), phi = c(1,1), sampleid = sampleid,spatial = T,
                                       intercept_unpenalized = TRUE,constr = F,
                                       maxit =100, tol = 1e-6, verbose = F)
      
      
      mu_full <- as.numeric(s_j * exp(X_full %*% fit_full$beta + ZS %*%fit_full$u))
      U = fit_full$u
      
      
      RobustTau = RobustTauPhi(as.vector(ZS %*%fit_full$u),mu_full,alpha,sampleid = sampleid)
      tau2 = pmin(pmax(1e-7,RobustTau$tau),10)
      phi  = pmin(pmax(1e-7,RobustTau$phi),10)
      ta   = ifelse(sum(tau2==0.01)>0, paste(which(tau2==0.01),collapse ="/"),0)
      
      fit_null <- fit_nb_glm_spatial_R(K = K, X = X_null, s = s_j, alpha = alpha,
                                       sigma2 = sigma2, coords = coords,Q_sparse=Q_sparse,spatialCov=spatialCov,
                                       tau2 = tau2, phi = phi, sampleid = sampleid,spatial = T,
                                       intercept_unpenalized = TRUE,constr =F,
                                       maxit = 100, tol = 1e-6, verbose = FALSE)
      
      
      mu_null <- as.numeric(s_j * exp(X_null %*% fit_null$beta+ ZS %*%fit_null$u))
      
      ll_null  <- nb_loglik(K , mu_null, alpha)
      
      
      
      fit_full <- fit_nb_glm_spatial_R(K = K, X = X_full, s = s_j, alpha = alpha,
                                       sigma2 = sigma2, coords = coords,Q_sparse=Q_sparse,spatialCov=spatialCov,
                                       tau2 = tau2, phi = phi, sampleid = sampleid,spatial = T,
                                       intercept_unpenalized = TRUE,constr =F,
                                       maxit =100, tol = 1e-6, verbose = F)
      
      mu_full <- as.numeric(s_j * exp(X_full %*% fit_full$beta +  ZS %*%fit_full$u))
      
      ll_full <- nb_loglik(K, mu_full, alpha)
      
      D <- 2 * (ll_full - ll_null)
      pval <- pchisq(D, df = ncol(X_full) - ncol(X_null), lower.tail = FALSE)
      
      wald_pvalue = 2 * pnorm(-abs(fit_full$beta[2]/sqrt(fit_full$cov[2,2])))
      
      c(stat = D, p = pval,gene =rownames(counts[i,,drop=F]),ta=ta )
    }
    
    res_list <- parLapply(cl, 1:G, test_gene)
    
    stopCluster(cl)
    
    res_mat <- do.call(rbind, res_list)
    gene = res_mat[, "gene"]
    res_mat <- lapply(1:2,function(i) res_mat[,i]%>%as.numeric) 
    res_mat <- do.call(cbind, res_mat)
    
    pvals <- res_mat[, 2]
    p.adj <- p.adjust(pvals, method = "BH")
    
    res_df <- data.frame(gene = gene, stat = res_mat[,1], pvalue = pvals, padj = p.adj)
    
    return(res_df)
    
  }
  
}



ClusterToTree <- function(Centers,weighted = TRUE){
  # Derive a tree from a center
  # INPUT: Centers
  # OUTPUT: mst (igraph) object.
  adjacency <- dist(Centers, method = "manhattan")
  full_graph <- graph.adjacency(as.matrix(adjacency), mode = "undirected", 
                                weighted = weighted)
  mst <- minimum.spanning.tree(full_graph)
  
  return(mst)
}


cauchy_combine <- function(pvals, weights = NULL){
  
  pvals <- pmin(pmax(pvals, 1e-15), 1-1e-15)
  
  if(is.null(weights)){
    weights <- rep(1/length(pvals), length(pvals))
  }
  
  Tstat <- sum(weights * tan((0.5 - pvals) * pi))
  pval  <- 0.5 - atan(Tstat) / pi
  
  return(pval)
}


computeSizeFactors <- function(countData) {
  # Step 1: compute geometric mean per gene
  geom_means <- apply(countData, 1, function(gene_counts) {
    # Avoid zeros by returning NA if all counts are zero
    if (all(gene_counts == 0)) return(NA)
    exp(mean(log(gene_counts[gene_counts > 0])))
  })
  
  # Step 2: compute ratios K_ij / g_i for each gene
  ratios <- countData / geom_means
  
  # Step 3: take median of ratios for each sample
  size_factors <- apply(ratios, 2, function(x) median(x[!is.na(x) & is.finite(x)], na.rm = TRUE))+0.0001
  
  # Normalize size factors to have geometric mean = 1 (optional)
  size_factors <- size_factors / exp(mean(log(size_factors+0.0001)))
  
  return(size_factors)
}


estimate_alpha_mom_vec <- function(counts, size_factors = NULL, min_mu = 1e-8, floor_alpha = 1e-8) {
  # counts: G x n matrix (genes x samples)
  if (is.null(size_factors)) size_factors <- rep(1, ncol(counts))
  G <- nrow(counts)
  alphas <- numeric(G)
  mu_bar <- numeric(G)
  
  for (i in seq_len(G)) {
    k <- as.numeric(counts[i, ])
    y <- k / size_factors
    m <- mean(y)
    v <- var(y)
    # Var(y) = m + alpha * m^2 → alpha = (v - m) / m^2
    if (m < min_mu) {
      alphas[i] <- floor_alpha
      mu_bar[i] <- m
    } else {
      a <- (v - m) / (m^2)
      if (!is.finite(a) || a < floor_alpha) a <- floor_alpha
      alphas[i] <- a
      mu_bar[i] <- m
    }
  }
  
  
  gene_names <- rownames(counts)
  if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(G))
  names(alphas) <- gene_names
  names(mu_bar) <- gene_names
  
  list(alpha = alphas, mu = mu_bar)
}


fit_dispersion_trend <- function(mu, alpha, use_loess = FALSE, loess_span = 0.5) {
  # mu, alpha are named numeric vectors
  ok <- is.finite(mu) & is.finite(alpha) & (mu > 0)
  mu0 <- mu[ok]; alpha0 <- alpha[ok]
  if (length(mu0) < 10) stop("Not enough genes to fit trend")
  
  if (use_loess) {
    lo <- loess(log(alpha0) ~ log(mu0), span = loess_span)
    alpha_trend_fn <- function(muq) exp(predict(lo, newdata = data.frame(mu0 = muq), se = FALSE))
    fitted <- alpha_trend_fn(mu)
    names(fitted) <- names(alpha)
    return(list(trend_fn = alpha_trend_fn,
                fitted = fitted,
                method = "loess"))
  } else {
    start_a0 <- min(median(alpha0), 1e-3)
    start_a1 <- max(median((alpha0 - start_a0) * mu0, na.rm = TRUE), 1e-6)
    df <- data.frame(mu = mu0, alpha = alpha0)
    
    fit_ok <- FALSE
    fit_par <- NULL
    try({
      nls_fit <- nls(alpha ~ a1/mu + a0, data = df,
                     start = list(a1 = start_a1, a0 = start_a0),
                     control = nls.control(maxiter = 200, warnOnly = TRUE))
      coefs <- coef(nls_fit)
      fit_ok <- TRUE
      fit_par <- coefs
    }, silent = TRUE)
    
    if (!fit_ok) {
      invmu <- 1 / mu0
      w <- 1 / pmax(alpha0, 1e-6)
      lm_fit <- lm(alpha0 ~ invmu, weights = w)
      coefs_lm <- coef(lm_fit)
      fit_par <- c(a0 = coefs_lm[1], a1 = coefs_lm[2])
    }
    
    alpha_trend_fn <- function(muq) {
      a1 <- fit_par["a1"]; a0 <- fit_par["a0"]
      a1 / muq + a0
    }
    fitted <- alpha_trend_fn(mu)
    names(fitted) <- names(alpha)
    
    return(list(trend_fn = alpha_trend_fn,
                fitted = fitted,
                params = fit_par,
                method = "parametric"))
  }
}

shrink_dispersions <- function(mu, alpha_gene, alpha_trend,
                               n_samples, n_params = 2,
                               sampling_var_method = "analytic", bootstrap_nboot = 0, s = NULL) {
  G <- length(alpha_gene)
  log_alpha_gene <- log(alpha_gene)
  log_alpha_trend <- log(alpha_trend)
  
  res_log <- log_alpha_gene - log_alpha_trend
  
  # estimate sampling variances for each gene: either analytic constant or bootstrap per gene
  if (sampling_var_method == "analytic") {
    sigma_i2 <- rep(approx_logalpha_sampling_var(n_samples, n_params, method = "analytic"), G)
  } else if (sampling_var_method == "constant") {
    sigma_i2 <- rep(bootstrap_nboot, G)
  } else {
    stop("Unknown sampling_var_method")
  }
  
  obs_var <- var(res_log, na.rm = TRUE)
  mean_sampling_var <- mean(sigma_i2, na.rm = TRUE)
  sigma_d2 <- obs_var - mean_sampling_var
  if (!is.finite(sigma_d2) || sigma_d2 < 0) sigma_d2 <- 0
  
  #Ensure vector length always matches number of genes
  weights <- rep(0, G)
  if (sigma_d2 > 0) {
    weights <- sigma_d2 / (sigma_d2 + sigma_i2)
  }
  
  log_alpha_final <- weights * log_alpha_gene + (1 - weights) * log_alpha_trend
  alpha_final <- exp(log_alpha_final)
  
  # Assign names consistently
  names(alpha_final) <- names(alpha_gene)
  names(weights) <- names(alpha_gene)
  names(res_log) <- names(alpha_gene)
  
  list(
    alpha_final = alpha_final,
    alpha_trend = alpha_trend,
    sigma_d2 = sigma_d2,
    weights = weights,
    sigma_i2 = sigma_i2,
    log_alpha_gene = log_alpha_gene,
    log_alpha_trend = log_alpha_trend,
    res_log = res_log
  )
}

approx_logalpha_sampling_var <- function(n_samples, n_params = 2, method = c("analytic","constant"), constant = NULL) {
  method <- match.arg(method)
  if (method == "analytic") {
    df <- n_samples - n_params
    if (df <= 0) stop("Not enough degrees of freedom for analytic approx")
    2 / df
  } else {
    if (is.null(constant)) stop("Provide constant for sampling var")
    constant
  }
}

fit_nb_glm_spatial_R <- function(K, X, s, alpha, sigma2,
                                 coords, Q_sparse,ZS, tau2 = 1, phi = 1,sampleid,spatial=T,
                                 intercept_unpenalized = TRUE,constr = T,
                                 maxit = 100, tol = 1e-8, verbose = FALSE) {
  n <- ncol(ZS )
  p <- ncol(X)
  
  # Fixed-effect prior
  Lambda_beta <- diag(1 / sigma2, p)
  if (intercept_unpenalized) Lambda_beta[1, 1] <- 0
  
  # Spatial covariance (Gaussian kernel)
  if(spatial){
    
    
    #### EXPONENTIAL
    #  uq = unique(sampleid)
    # QQ = list()
    # for (l in 1:length(uq)) {
    
    #  D <- as.matrix(dist(coords[sampleid==uq[l],]))
    #  Sigma <- exp(- (D^2) / (2 * phi[l]^2))
    
    
    ## MATERN COVRIANCE
    # nugget_eps <- 1e-6 
    # uq = unique(sampleid)
    # nu=rep(1.5,length(uq))
    # QQ = list()
    #  for (l in 1:length(uq)) {
    
    #    D <- as.matrix(dist(coords[sampleid==uq[l],]))
    #rho_l <- phi[l] 
    #Sigma <- fields::Matern(D, range = rho_l, smoothness = nu[l])
    #diag(Sigma) <- diag(Sigma) + nugget_eps
    #chol_S <- tryCatch(chol(Sigma), error = function(e) NULL)
    #if (!is.null(chol_S)) {
    #  Sigma_inv <- chol2inv(chol_S)
    #} else {
    # fallback to a more stable solver (regularize a bit)
    #  Sigma_reg <- Sigma + diag(nugget_eps * 10, nrow(Sigma))
    #   Sigma_inv <- tryCatch(solve(Sigma_reg), error = function(e) {
    # last fallback: pseudo-inverse
    #   as.matrix(MASS::ginv(Sigma_reg))
    #  })
    #}
    #QQ[[l]] <- Sigma_inv / tau2[l]
    
    
    
    
    # Prior precision for spatial random effects
    #   QQ[[l]] <- solve(Sigma + diag(1e-6, nrow(Sigma))) / tau2[l]
    
    # }
    
    
    
    if(constr & ncol(X)>1){ 
      spatialCov = as.vector(ZS%*%matrix(1:ncol(ZS)))
      Xconsr = data.frame(X[,-1,drop=F],scov = spatialCov) %>%
        group_by(scov)%>%  dplyr::summarise_all(mean) %>% 
        dplyr::select(-scov) %>%as.matrix()
      
      P_X <- Xconsr %*% solve(t(Xconsr) %*% Xconsr) %*% t(Xconsr)
      H <- diag(nrow(Xconsr)) - P_X
      
      tmp <- Q_sparse %*% H
      Q_sparse <- crossprod(H, tmp)
      #Q_sparse  <- H %*% Q_sparse %*% H 
    }
    
    X_star <- cbind(X, ZS)
    
    
    Lambda_star <- bdiag(Lambda_beta, Q_sparse)
    
    # Initial parameters
    beta_u <- rep(0, p + n)
    
    for (iter in 1:maxit) {
      # eta <- as.vector(log(s) + X_star %*% beta_u)
      # mu <- exp(eta)
      eta <- as.vector(log(s + 1e-8) + X_star %*% beta_u)
      mu <- exp(pmin(eta, 14))  # cap log-scale at exp(14) ~ 1.2e6
      
      # NB variance: var = mu + alpha * mu^2
      W <- 1 / (1/mu + alpha)
      z <- eta + (K - mu) / mu
      z[!is.finite(z)] <- eta[!is.finite(z)]
      z <- pmax(pmin(z, 20), -20)  # keep pseudo-response within ±20
      
      WX <- X_star * sqrt(W)
      wz <- sqrt(W) * z
      
      ridge <- 1e-5 * diag(ncol(WX))
      
      WX <- Matrix(WX, sparse = TRUE)
      H <- crossprod(WX) + as.matrix(Lambda_star)+ridge
      rhs <- crossprod(WX, wz)
      
      #beta_new <- solve(H, rhs)
      
      beta_candidate <- solve(H, rhs)
      delta <- beta_candidate - beta_u
      
      if (any(abs(delta) > 10)) {  # limit huge jumps
        delta <- sign(delta) * pmin(abs(delta), 10)
      }
      beta_new <- beta_u + delta
      
      beta_new[(p + 1):(p + n)] = beta_new[(p + 1):(p + n)]-mean(beta_new[(p + 1):(p + n)])
      #beta_new[(p + 1):(p + n/2)] = beta_new[(p + 1):(p + n/2)]-mean(beta_new[(p + 1):(p + n/2)])
      #beta_new[((p + 1):(p + n/2))+n/2 ] = beta_new[((p + 1):(p + n/2))+n/2 ]-mean(beta_new[((p + 1):(p + n/2))+n/2 ])
      
      diff <- max(abs(beta_new - beta_u),na.rm = T)
      if (verbose) message("Iter ", iter, " | Δ=", round(diff, 6))
      
      if(is.nan(diff)|is.na(diff)) break
      
      beta_u <- beta_new
      
      if(is.nan(diff)|is.na(diff)){
        warning("Input is empty, returning NULL.", call. = FALSE)
        return(NULL)
      }
      
      if (diff < tol) break
    }
    
  }else{
    
    # Augmented design matrix
    X_star <- cbind(X)
    Lambda_star <- bdiag(Lambda_beta)
    
    # Initial parameters
    beta_u <- rep(0, p)
    
    for (iter in 1:maxit) {
      # eta <- as.vector(log(s) + X_star %*% beta_u)
      # mu <- exp(eta)
      eta <- as.vector(log(s + 1e-8) + X_star %*% beta_u)
      mu <- exp(pmin(eta, 14))  # cap log-scale at exp(14) ~ 1.2e6
      
      # NB variance: var = mu + alpha * mu^2
      W <- 1 / (1/mu + alpha)
      z <- eta + (K - mu) / mu
      z[!is.finite(z)] <- eta[!is.finite(z)]
      z <- pmax(pmin(z, 20), -20)  # keep pseudo-response within ±20
      
      WX <- X_star * sqrt(W)
      wz <- sqrt(W) * z
      
      ridge <- 1e-8 * diag(ncol(WX))
      
      H <- crossprod(WX) + as.matrix(Lambda_star)+ridge
      rhs <- crossprod(WX, wz)
      
      #beta_new <- solve(H, rhs)
      
      beta_candidate <- solve(H, rhs)
      delta <- beta_candidate - beta_u
      
      if (any(abs(delta) > 10)) {  # limit huge jumps
        delta <- sign(delta) * pmin(abs(delta), 10)
      }
      beta_new <- beta_u + delta
      
      diff <- max(abs(beta_new - beta_u),na.rm = T)
      if (verbose) message("Iter ", iter, " | Δ=", round(diff, 6))
      
      if(is.nan(diff)|is.na(diff)) break
      
      beta_u <- beta_new
      if(is.nan(diff)|is.na(diff)){
        warning("Input is empty, returning NULL.", call. = FALSE)
        return(NULL)
      }
      if (diff < tol) break
      
    }
    
  }
  
  
  cov_beta_u <- solve(crossprod(WX) + as.matrix(Lambda_star)+ridge)
  se_beta_u <- sqrt(diag(cov_beta_u))
  
  list(beta = beta_u[1:p],
       u = beta_u[(p + 1):(p + n)],
       se_beta = se_beta_u[1:p],
       se_u = se_beta_u[(p + 1):(p + n)],
       cov = cov_beta_u,
       iterations = iter,
       converged = diff < tol)
}


RobustTauPhi <- function(K,mu_full,alpha,sampleid){
  
  #r = (K-mu_full)/(sqrt(mu_full+mu_full^2*alpha)+0.0001)
  r = (K)/(sqrt(mu_full+mu_full^2*alpha)+0.0001)
  uq = unique(sampleid)
  d=NULL
  for (i in uq) {
    id = sampleid==i
    
    df <- data.frame(x = coords[id,1], y = coords[id,2], res = r[id])
    df_sf <- st_as_sf(df, coords = c("x", "y"), crs = NA)
    D <- as.matrix(dist(coords))
    bin_width <- max(D) / 10
    vg <- variogram(res ~ 1, data = df_sf, width = bin_width)
    model_fun <- function(par, d) {
      nug <- par[1]; sill <- par[2]; phi <- par[3]
      nug + sill * (1 - exp(-d / phi))
    }
    loss_fun <- function(par) {
      pred <- model_fun(par, vg$dist)
      w <- sqrt(vg$np)  # weight by sqrt(np) or np
      sum(w * (vg$gamma - pred)^2)
    }
    init <- c(nugget=0, sill=var(r)*0.5, phi=median(vg$dist))
    opt <- optim(init, loss_fun, method="L-BFGS-B", lower=c(0,0, min(vg$dist)/10), upper=c(Inf,Inf, max(vg$dist)*10))
    par_hat <- opt$par
    
    ds = data.frame(
      tau = opt$par[2],
      phi = opt$par[3],
      nugget=opt$par[1]
    )
    rownames(ds) = paste0("Grp",i)
    d = rbind(d,ds)
  }
  return(d)
}


RobustTauPhi_Matern <- function(K, mu_full, alpha, sampleid, coords, nu = 1.5) {
  # Residuals standardized by NB variance
  r <- K#(K - mu_full) / (sqrt(mu_full + mu_full^2 * alpha) + 1e-4)
  r <- K/ (sqrt(mu_full + mu_full^2 * alpha) + 1e-4)
  
  uq <- unique(sampleid)
  d <- NULL
  
  for (i in uq) {
    id <- sampleid == i
    
    # Spatial coordinates + residuals for this group
    df <- data.frame(x = coords[id, 1], y = coords[id, 2], res = r[id])
    df_sf <- st_as_sf(df, coords = c("x", "y"), crs = NA)
    
    # Compute empirical variogram
    D <- as.matrix(dist(coords[id, ]))
    bin_width <- max(D) / 10
    vg <- variogram(res ~ 1, data = df_sf, width = bin_width)
    
    # --- Matérn model function (correlation, normalized) ---
    model_fun <- function(par, d) {
      nug <- par[1]
      sill <- par[2]
      phi  <- par[3]
      # Matérn covariance with unit variance and smoothness nu
      # fields::Matern() can compute this directly; here we write explicitly:
      arg <- d / phi
      matern_corr <- (2^(1 - nu)) / gamma(nu) * (arg^nu) * besselK(arg, nu)
      matern_corr[arg == 0] <- 1  # handle diagonal
      # Convert covariance to semivariogram
      nug + sill * (1 - matern_corr)
    }
    
    # Weighted least squares loss (robust)
    loss_fun <- function(par) {
      pred <- model_fun(par, vg$dist)
      w <- sqrt(vg$np)
      sum(w * (vg$gamma - pred)^2, na.rm = TRUE)
    }
    
    # Initialize parameters
    init <- c(nugget = 0.01 * var(r[id], na.rm = TRUE),
              sill   = 0.5 * var(r[id], na.rm = TRUE),
              phi    = median(vg$dist, na.rm = TRUE))
    
    # Optimize
    opt <- optim(
      par = init,
      fn = loss_fun,
      method = "L-BFGS-B",
      lower = c(0, 0, min(vg$dist) / 10),
      upper = c(Inf, Inf, max(vg$dist) * 10)
    )
    
    # Store result
    ds <- data.frame(
      tau = opt$par[2],        # sill (spatial variance)
      phi = opt$par[3],        # range
      nugget = opt$par[1],
      group = i
    )
    rownames(ds) <- paste0("Grp", i)
    d <- rbind(d, ds)
  }
  return(d)
}



fit_binorm_glm_spatial_R <- function(K, X,sigma2,
                                     coords,Q_sparse,ZS, tau2 = 1, phi = 1,sampleid,spatial=T,
                                     intercept_unpenalized = TRUE,
                                     maxit = 100, tol = 1e-8, verbose = FALSE) {
  n <- ncol(ZS)
  p <- ncol(X)
  
  # Fixed-effect prior
  Lambda_beta <- diag(1 / sigma2, p)
  if (intercept_unpenalized) Lambda_beta[1, 1] <- 0
  
  # Spatial covariance (Gaussian kernel)
  if(spatial){
    
    
    # nugget_eps <- 1e-6 
    # uq = unique(sampleid)
    # nu = rep(1.5,length(uq))
    #QQ = list()
    # for (l in 1:length(uq)) {
    
    #   D <- as.matrix(dist(coords[sampleid==uq[l],]))
    
    #   rho_l <- phi[l] 
    
    ## MATERN COVRIANCE
    
    #  Sigma <- fields::Matern(D, range = rho_l, smoothness = nu[l])
    #  diag(Sigma) <- diag(Sigma) + nugget_eps
    #  chol_S <- tryCatch(chol(Sigma), error = function(e) NULL)
    #  if (!is.null(chol_S)) {
    #    Sigma_inv <- chol2inv(chol_S)
    #  } else {
    # fallback to a more stable solver (regularize a bit)
    #    Sigma_reg <- Sigma + diag(nugget_eps * 10, nrow(Sigma))
    #   Sigma_inv <- tryCatch(solve(Sigma_reg), error = function(e) {
    #      # last fallback: pseudo-inverse
    #      as.matrix(MASS::ginv(Sigma_reg))
    #    })
    #  }
    #  QQ[[l]] <- Sigma_inv / tau2[l]
    
    # EXPONENTIAL
    
    #Sigma <- exp(- (D^2) / (2 * phi[l]^2))
    
    # Prior precision for spatial random effects
    #QQ[[l]] <- solve(Sigma + diag(1e-6, nrow(Sigma))) / tau2[l]
    
    # }
    #Q = do.call(bdiag,QQ)
    
    # Augmented design matrix
    #ZS = as.matrix(fastDummies::dummy_cols(spatialCov,remove_selected_columns = TRUE))
    
    X_star <- cbind(X, ZS)
    Lambda_star <- bdiag(Lambda_beta, Q_sparse)
    
    # Initial parameters
    beta_u <- rep(0.001, p + n)
    
    for (iter in 1:maxit) {
      # eta <- as.vector(log(s) + X_star %*% beta_u)
      # mu <- exp(eta)
      eta <- as.vector(X_star %*% beta_u)
      mu <- 1/(1+exp(-pmin(eta, 14)))  # cap log-scale at exp(14) ~ 1.2e6
      
      # NB variance: var = mu + alpha * mu^2
      W <- mu*(1-mu)
      z <- eta + (K - mu) / W
      z[!is.finite(z)] <- eta[!is.finite(z)]
      z <- pmax(pmin(z, 20), -20)  # keep pseudo-response within ±20
      
      WX <- X_star * sqrt(W)
      wz <- sqrt(W) * z
      
      ridge <- 1e-3 * diag(ncol(WX))
      
      WX <- Matrix(WX, sparse = TRUE)
      H <- crossprod(WX) + as.matrix(Lambda_star)+ridge
      rhs <- crossprod(WX, wz)
      
      #beta_new <- solve(H, rhs)
      
      beta_candidate <- solve(H, rhs)
      delta <- beta_candidate - beta_u
      
      if (any(abs(delta) > 10)) {  # limit huge jumps
        delta <- sign(delta) * pmin(abs(delta), 10)
      }
      beta_new <- beta_u + delta
      #beta_new[(p + 1):(p + n)] = beta_new[(p + 1):(p + n)]-mean(beta_new[(p + 1):(p + n)])
      beta_new[(p + 1):(p + n/2)] = beta_new[(p + 1):(p + n/2)]-mean(beta_new[(p + 1):(p + n/2)])
      beta_new[((p + 1):(p + n/2))+n/2 ] = beta_new[((p + 1):(p + n/2))+n/2 ]-mean(beta_new[((p + 1):(p + n/2))+n/2 ])
      diff <- max(abs(beta_new - beta_u),na.rm = T)
      if (verbose) message("Iter ", iter, " | Δ=", round(diff, 6))
      
      if(is.nan(diff)|is.na(diff)) break
      
      beta_u <- beta_new
      
      if(is.nan(diff)|is.na(diff)){
        warning("Input is empty, returning NULL.", call. = FALSE)
        return(NULL)
      }
      
      if (diff < tol) break
    }
    
  }else{
    
    # Augmented design matrix
    X_star <- cbind(X)
    Lambda_star <- bdiag(Lambda_beta)
    
    # Initial parameters
    beta_u <- rep(0, p)
    
    for (iter in 1:maxit) {
      # eta <- as.vector(log(s) + X_star %*% beta_u)
      # mu <- exp(eta)
      
      eta <- as.vector(X_star %*% beta_u)
      mu <- 1/(1+exp(-pmin(eta, 14)))  # cap log-scale at exp(14) ~ 1.2e6
      
      # NB variance: var = mu + alpha * mu^2
      W <- mu*(1-mu)
      z <- eta + (K - mu) / W
      z[!is.finite(z)] <- eta[!is.finite(z)]
      z <- pmax(pmin(z, 20), -20)  # keep pseudo-response within ±20
      
      WX <- X_star * sqrt(W)
      wz <- sqrt(W) * z
      
      ridge <- 1e-4 * diag(ncol(WX))
      
      H <- crossprod(WX) + as.matrix(Lambda_star)+ridge
      rhs <- crossprod(WX, wz)
      
      #beta_new <- solve(H, rhs)
      
      beta_candidate <- solve(H, rhs)
      delta <- beta_candidate - beta_u
      
      if (any(abs(delta) > 10)) {  # limit huge jumps
        delta <- sign(delta) * pmin(abs(delta), 10)
      }
      beta_new <- beta_u + delta
      
      diff <- max(abs(beta_new - beta_u),na.rm = T)
      if (verbose) message("Iter ", iter, " | Δ=", round(diff, 6))
      
      if(is.nan(diff)|is.na(diff)) break
      
      beta_u <- beta_new
      if(is.nan(diff)|is.na(diff)){
        warning("Input is empty, returning NULL.", call. = FALSE)
        return(NULL)
      }
      if (diff < tol) break
      
    }
    
  }
  
  
  cov_beta_u <- solve(crossprod(WX) + as.matrix(Lambda_star)+ridge)
  se_beta_u <- sqrt(diag(cov_beta_u))
  
  list(beta = beta_u[1:p],
       u = beta_u[(p + 1):(p + n)],
       se_beta = se_beta_u[1:p],
       se_u = se_beta_u[(p + 1):(p + n)],
       cov = cov_beta_u,
       iterations = iter,
       converged = diff < tol)
}


book.mesh.dual <- function(mesh) {
  # Function to construct dual Mesh
  
  #Input: mesh (inla mesh object)
  #Output: dual mesh (inla object)
  # Referece -> INLA package
  
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

