#  Spatial-ZEDNet: zero-inflated graphical modeling of exposure-induced differential expression in spatial transcriptomics
Differential expression analysis in spatially resolved single-cell data enables the detection of exposure-induced spatial patterns and localized gene regulation. While numerous methods exist for identifying differentially expressed genes (DEGs) in conventional single-cell RNA sequencing (scRNA-seq), few explicitly incorporate spatial context. A major challenge in spatial DEG analysis is the misalignment of tissue sections between control and perturbed conditions, making coordinate-based comparisons unreliable. Furthermore, single-cell data often exhibit zero inflation, a high frequency of zero counts due to technical dropouts, which can obscure true biological differences when not properly modeled. To address these challenges, we present Spatial-ZEDNet, a novel framework for identifying exposure-induced DEGs in spatially structured single-cell data from two conditions. Spatial-ZEDNet integrates spatial information through a hierarchical network-based Gaussian field model that accounts for spatial variability, models zero-inflated gene expression distributions, and leverages tissue architecture to align biological signals across exposed and control conditions. We validated our method through simulation studies and benchmarking against existing tools. Application to two real datasets, colitis-induced tissue remodeling, and Plasmodium infection, revealed spatially organized gene expression patterns and key pathways implicated in disease progression. Our findings demonstrate that Spatial-ZEDNet is a powerful and interpretable tool for spatial transcriptomic analysis and may facilitate the development of targeted therapies for inflammatory and infectious diseases.

![image](https://github.com/user-attachments/assets/36cb2e1b-41ed-4426-ba43-c85d3b8212d4)

# Demo analysis

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
                          countmodel = "lognormal", # nbinomial for count,
                          CollectPostD = TRUE)
```

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

