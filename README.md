#  Spatial-ZEDNet: zero-inflated graphical modeling of exposure-induced differential expression in spatial transcriptomics
Differential expression analysis in spatially resolved single-cell data enables the detection of exposure-induced spatial patterns and localized gene regulation. While numerous methods exist for identifying differentially expressed genes (DEGs) in conventional single-cell RNA sequencing (scRNA-seq), few explicitly incorporate spatial context. A major challenge in spatial DEG analysis is the misalignment of tissue sections between control and perturbed conditions, making coordinate-based comparisons unreliable. Furthermore, single-cell data often exhibit zero inflation, a high frequency of zero counts due to technical dropouts, which can obscure true biological differences when not properly modeled. To address these challenges, we present Spatial-ZEDNet, a novel framework for identifying exposure-induced DEGs in spatially structured single-cell data from two conditions. Spatial-ZEDNet integrates spatial information through a hierarchical network-based Gaussian field model that accounts for spatial variability, models zero-inflated gene expression distributions, and leverages tissue architecture to align biological signals across exposed and control conditions. We validated our method through simulation studies and benchmarking against existing tools. Application to two real datasets, colitis-induced tissue remodeling, and Plasmodium infection, revealed spatially organized gene expression patterns and key pathways implicated in disease progression. Our findings demonstrate that Spatial-ZEDNet is a powerful and interpretable tool for spatial transcriptomic analysis and may facilitate the development of targeted therapies for inflammatory and infectious diseases.

![image](https://github.com/user-attachments/assets/36cb2e1b-41ed-4426-ba43-c85d3b8212d4)

# Demo analysis
```{R}
# Required Packages
library(INLA) # For Imputation
library(tidyverse) # Data manipulation
library(igraph)
library(scales)
library(genie)
library(ggraph)
library(progress)
library(viridis)

library(parallel)
library(parallelly)
library(progress)
library(progressr)
library(sp)
library(sf)
library(MASS)
library(gstat)

```
## Find differential and activated genes

```{R}
# Load data
library(Seurat)
SeuratObj = readRDS("SampleData.rds")

# Get utility function

source("https://raw.githubusercontent.com/NIEHS/Spatial-ZEDINet/main/SpatialZEDNet_Git.R")

# Find differential and activated Genes

AuxResult =  SpatialDENet(all_data = SeuratObj@assays$RNA$data,
                          metadata = SeuratObj@meta.data, # Must contain "Group" column
                          Sample_id = "Sample_id",
                          joint = TRUE,
                          countmodel = "lognormal", # nbinomial for count,
                          CollectPostD = TRUE)
AuxResult$adj_pvalue_BY
```
![image](https://github.com/user-attachments/assets/2c365625-cc78-47ea-85e6-2c8038c05ba8)
## Plots 
```{R}
# Plot cell types on the network

Metadata = AuxResult$Metadata
Metadata$cellType =  as.numeric(as.factor(as.character(Metadata$Tier1 )))
Pl = Metadata %>%group_by(meshid_id) %>%summarise(celltype=round(mean(cellType)))
Pl$celltype_lab = factor(Pl$celltype,levels = 1:4,labels = c("Endothelial","Epithelial","Fibroblast","Immune"))

plotTree(mst,Pl$celltype_lab,cols = c25[8:12],Lab = F,
         edge_alpha = .01,vertex.size =AuxResult$Res$nn,
         legend.size = 3)
```

![image](https://github.com/user-attachments/assets/b3de9ba1-8ca5-4537-8c4e-3b4128e37d7b)



```{R}
#### Disease and health

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


Pl = Metadata %>%group_by(meshid_id) %>%summarise(type=Mode(Sample_type))


plotTree(mst,Pl$type,cols = c25[3:4],Lab = F,
         edge_alpha = .01,vertex.size =AuxResult$Res$nn,
         legend.size = 3)
```

![image](https://github.com/user-attachments/assets/150ab57e-bc95-47b4-9624-ce5d6b0c477a)


```{R}
adj_pvalue = AuxResult$adj_pvalue_BY[1,]
nam = names(sort(adj_pvalue[adj_pvalue <0.01],decreasing = F))
mst = AuxResult$mst

##### Network effect Count
antibody= "Lag3"
plotTree(mst,AuxResult$treeEffect_cont[,antibody],
         edge_alpha = .01,
         Lab = F,
         vertex.size =AuxResult$Res$nn,
         main=antibody)
```
![image](https://github.com/user-attachments/assets/606f622d-259c-47ae-9f66-f004650feaac)
```{R}
##### Network effect Activation
plotTree(mst,AuxResult$treeEffect_actv[,antibody],
         edge_alpha = .01,
         Lab = F,
         vertex.size =AuxResult$Res$nn,
         main= antibody)

```
![image](https://github.com/user-attachments/assets/c047ae6b-9f98-4ac4-911b-404067e32d92)

```{R}
### Plot effect on image

Centers = data.frame(meshid_id = 1:vcount(mst), AuxResult$treeEffect_cont)
Metadata =AuxResult$Metadata
Metadata = left_join(Metadata,Centers,by= "meshid_id")

mn = min(Metadata[,antibody])
ma = max(Metadata[,antibody])

Pl = Metadata %>% filter(Slice_ID =="062921_D0_m3a_2_slice_3")

p1 = plotScatter(Pl$x,Pl$y,Pl[,antibody],main = paste0("Day 0: ",antibody), ManualColor = F,cols = c25[8:12],
            legend.size = 3,limits = c(mn,ma))


Pl = Metadata %>% filter(Slice_ID =="062221_D9_m3_2_slice_3")

p2= plotScatter(Pl$x,Pl$y,Pl[,antibody],main = paste0("Day 9: ",antibody),ManualColor = F,cols = c25[8:12],
            legend.size = 3,limits = c(mn,ma))

p1|p2
```

![image](https://github.com/user-attachments/assets/b643659b-4c67-430c-91e4-3fb7cb018a28)

```{R}
# Plot Corregulation matrix
nam_Dnet_nam = nam

M = Metadata[,c("Sample_type",nam_Dnet_nam[(nam_Dnet_nam%in%colnames(Metadata))])]
M = M %>%filter(Sample_type=="DSS9") %>%dplyr::select(-Sample_type)
M = unique(M) %>%cor

pheatmap::pheatmap(M, 
                   treeheight_row = 0,
                   treeheight_col = 0,
                   fontsize = 4.5,
                   angle_col = 90,
                   cluster_rows = T,
                   cluster_cols = T,
                   #breaks = breaks,
                   color = inferno(49))
```
![image](https://github.com/user-attachments/assets/ae29a13b-f3fa-43ca-958d-fdd8513cd0af)

```{R}
####### Plot effect on matrix 
M = Metadata[,c("Sample_type",nam_Dnet_nam[(nam_Dnet_nam%in%colnames(Metadata))])]
M = M %>%filter(Sample_type=="Healthy") %>%dplyr::select(-Sample_type)
M1 = unique(M)

M = Metadata[,c("Sample_type",nam_Dnet_nam[(nam_Dnet_nam%in%colnames(Metadata))])]
M = M %>%filter(Sample_type=="DSS9") %>%dplyr::select(-Sample_type)
M2 = M%>%unique()

# Cluster M1 and M2 separately
hc_M1 <- hclust(dist(M1))
hc_M2 <- hclust(dist(M2))
M1_clustered <- M1[hc_M1$order, ]
M2_clustered <- M2[hc_M2$order, ]
combined_mat <- rbind(M1_clustered, M2_clustered)

# Row annotation
annotation_row <- data.frame(Source = c(rep("Day 0", nrow(M1)), rep("Day 9", nrow(M2))))
rownames(annotation_row) <- rownames(combined_mat)
purple_yellow_palette <- colorRampPalette(c("#762A83", "#FFFFFF", "#d95f02"))(100)
ann_colors <- list(Source = c(`Day 0` = "#1b9e77", `Day 9` = "#FFDC00")) 
# Heatmap
pheatmap(combined_mat,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         scale = "row",
         color = purple_yellow_palette,
         fontsize = 5,
         treeheight_col = 0
)
```

![image](https://github.com/user-attachments/assets/2e5a87f1-7a9f-4ebd-84ed-52096a88f5de)

