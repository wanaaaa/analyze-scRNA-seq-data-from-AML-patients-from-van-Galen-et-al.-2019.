library(readr)
library(tidyverse)
library(Seurat)
library(ggsci)
library(readr)
library(dplyr)

SeuratFun <- function(dem_tbl, annot_tbl) {
  print('In SeuratFun')
  dem_mat <- dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
  dem_mat[1:5, 1:3]
  # Check the number of UMIs/counts per cell/column.
  dem_mat %>% colSums() %>% summary()
  
  # The annotation data should be a data frame.
  anno_df <- annot_tbl %>% column_to_rownames("Cell") %>% as.data.frame()
  
  # ======================================================
  # Create the Seurat object
  # ======================================================
  
  # Initialize the Seurat object from the raw counts matrix.
  Seurat_so <- CreateSeuratObject(counts = dem_mat, min.cells = 1, min.features = 1)
  
  # The [[ operator can add columns to object metadata. We can add the percentage of reads that map to the mitochondrial genome.
  Seurat_so[["pct_mito"]] <- PercentageFeatureSet(Seurat_so, pattern = "^MT-")

  # add additional metadata using AddMetaData().
  Seurat_so <- AddMetaData(Seurat_so, metadata = anno_df)
  # head(Seurat_so@meta.data)

  # These cells are already pre-filtered, but we can visualize the QC metrics anyway.
  VlnPlot(Seurat_so, features = c("nFeature_RNA", "nCount_RNA", "pct_mito"), ncol = 3)

  # Perform standard analysis (log-transformed)

  # The default global-scaling normalization method LogNormalize normalizes the feature expression measurements for each cell
  # by the total expression, multiplies by a scale factor of 10,000, and log-transforms the result.
  Seurat_so <- NormalizeData(Seurat_so , normalization.method = "LogNormalize")
  
  return(Seurat_so)
}


grouped_soFun <- function(iSeurat_so) {
  # # # Identify highly variable features, scale, and perform dimensionality reduction.
  Seurat_so <- FindVariableFeatures(iSeurat_so, verbose = FALSE)
  Seurat_so <- ScaleData(Seurat_so, vars.to.regress = "pct_mito", verbose = FALSE)
  Seurat_so <- RunPCA(Seurat_so, verbose = FALSE)
  Seurat_so <- RunUMAP(Seurat_so, dims = 1:30, verbose = FALSE)
  
  return(Seurat_so)
}


BM3_dem_tbl <- read_tsv("../Data/hwData/GSM3587998_BM3.dem.txt.gz", show_col_types = FALSE)
BM3_annot_tbl <- read_tsv("../Data/hwData/GSM3587999_BM3.anno.txt.gz", show_col_types = FALSE)
BM3_so = SeuratFun(BM3_dem_tbl, BM3_annot_tbl)


BM4_dem_tbl <- read_tsv("../Data/hwData/GSM3588000_BM4.dem.txt.gz", show_col_types = FALSE)
BM4_annot_tbl <- read_tsv("../Data/hwData/GSM3588001_BM4.anno.txt.gz", show_col_types = FALSE)
BM4_so = SeuratFun(BM4_dem_tbl, BM4_annot_tbl)


BM5.34p_dem_tbl <- read_tsv("../Data/hwData/GSM3588002_BM5-34p.dem.txt.gz", show_col_types = FALSE)
BM5.34p_annot_tbl <- read_tsv("../Data/hwData/GSM3588002_BM5-34p.anno.txt.gz", show_col_types = FALSE)
BM5.34p_so = SeuratFun(BM5.34p_dem_tbl, BM5.34p_annot_tbl)


BM5.34p38n_dem_tbl <- read_tsv("../Data/hwData/GSM3588003_BM5-34p38n.dem.txt.gz", show_col_types = FALSE)
BM5.34p38n_annot_tbl <- read_tsv("../Data/hwData/GSM3588003_BM5-34p38n.anno.txt.gz", show_col_types = FALSE)
BM5.34p38n_so = SeuratFun(BM5.34p38n_dem_tbl, BM5.34p38n_annot_tbl)

AML419_dem_tbl <- read_tsv("../Data/hwData/GSM3587950_AML419A-D0.dem.txt.gz", show_col_types = FALSE)
AML419_annot_tbl <- read_tsv("../Data/hwData/GSM3587951_AML419A-D0.anno.txt.gz", show_col_types = FALSE)
AML419_so = SeuratFun(AML419_dem_tbl, AML419_annot_tbl)


AML707_dem_tbl <- read_tsv("../Data/hwData/GSM3587969_AML707B-D0.dem.txt.gz", show_col_types = FALSE)
AML707_annot_tbl <- read_tsv("../Data/hwData/GSM3587970_AML707B-D0.anno.txt.gz", show_col_types = FALSE)
AML707_so = SeuratFun(AML707_dem_tbl, AML707_annot_tbl)


AML916_dem_tbl <- read_tsv("../Data/hwData/GSM3587988_AML916-D0.dem.txt.gz", show_col_types = FALSE)
AML916_annot_tbl <- read_tsv("../Data/hwData/GSM3587989_AML916-D0.anno.txt.gz", show_col_types = FALSE)
AML916_so = SeuratFun(AML916_dem_tbl, AML916_annot_tbl)


AML921_dem_tbl <- read_tsv("../Data/hwData/GSM3587990_AML921A-D0.dem.txt.gz", show_col_types = FALSE)
AML921_annot_tbl <- read_tsv("../Data/hwData/GSM3587991_AML921A-D0.anno.txt.gz", show_col_types = FALSE)
AML921_so = SeuratFun(AML921_dem_tbl, AML921_annot_tbl)

# merge Seurat objects
BM_so = merge(x=BM3_so, y= c(BM4_so, BM5.34p_so, BM5.34p38n_so))
# merge Seurat objects
AML_so = merge(x=AML419_so, y= c(AML707_so, AML916_so, AML921_so))
# merge Seurat objects
BM.AML_so <- merge(BM_so, AML_so)

# grouping process with PCA, UMAP
groupedBM_so <- grouped_soFun(BM_so)
groupedAML_so <- grouped_soFun(AML_so)

DimPlot(groupedAML_so, reduction = "umap", group.by = "CellType", shuffle = TRUE) + 
  theme(aspect.ratio = 1) + scale_color_igv() + ggtitle("AML by CellType")
DimPlot(groupedAML_so, reduction = "umap", group.by = "PredictionRefined", shuffle = TRUE) + theme(aspect.ratio = 1) + 
  scale_color_igv() + ggtitle("AML by Malignancy status")
DimPlot(groupedAML_so, reduction = "umap", shuffle = TRUE) + theme(aspect.ratio = 1) + 
  scale_color_igv() + ggtitle(" AML By patinet ")
