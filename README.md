#  Spatial-ZEDNet: zero-inflated graphical modeling of exposure-induced differential expression in spatial transcriptomics
Differential expression analysis in spatially resolved single-cell data enables the detection of exposure-induced spatial patterns and localized gene regulation. While numerous methods exist for identifying differentially expressed genes (DEGs) in conventional single-cell RNA sequencing (scRNA-seq), few explicitly incorporate spatial context. A major challenge in spatial DEG analysis is the misalignment of tissue sections between control and perturbed conditions, making coordinate-based comparisons unreliable. Furthermore, single-cell data often exhibit zero inflation, a high frequency of zero counts due to technical dropouts, which can obscure true biological differences when not properly modeled. To address these challenges, we present Spatial-ZEDNet, a novel framework for identifying exposure-induced DEGs in spatially structured single-cell data from two conditions. Spatial-ZEDNet integrates spatial information through a hierarchical network-based Gaussian field model that accounts for spatial variability, models zero-inflated gene expression distributions, and leverages tissue architecture to align biological signals across exposed and control conditions. We validated our method through simulation studies and benchmarking against existing tools. Application to two real datasets, colitis-induced tissue remodeling, and Plasmodium infection, revealed spatially organized gene expression patterns and key pathways implicated in disease progression. Our findings demonstrate that Spatial-ZEDNet is a powerful and interpretable tool for spatial transcriptomic analysis and may facilitate the development of targeted therapies for inflammatory and infectious diseases.

![image](https://github.com/user-attachments/assets/36cb2e1b-41ed-4426-ba43-c85d3b8212d4)
```{R}
# Load data
library(Seurat)
SeuratObj = readRDS("SampleData.rds")
# Find differential and activated Genes

AuxResult =  SpatialDENet(all_data = SeuratObj@assays$RNA$data,
                          metadata = SeuratObj@meta.data, # Must contain "Group" column
                          Sample_id = "Sample_id",
                          countmodel = "lognormal", # nbinomial for count,
                          CollectPostD = TRUE)
```
