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
  # ======================================================
  
  # The default global-scaling normalization method LogNormalize normalizes the feature expression measurements for each cell
  # by the total expression, multiplies by a scale factor of 10,000, and log-transforms the result.
  Seurat_so <- NormalizeData(Seurat_so , normalization.method = "LogNormalize")
  
  # The normalized values are stored in the data slot in the RNA assay
  Seurat_so@assays$RNA$data[1:8, 1:3]
  # Each cell now has 10,000 normalized counts.
  Seurat_so@assays$RNA$data %>% expm1() %>% colSums() %>% summary()
  # 
  # # Identify highly variable features, scale, and perform dimensionality reduction.
  Seurat_so <- FindVariableFeatures(Seurat_so, verbose = FALSE)
  Seurat_so <- ScaleData(Seurat_so, vars.to.regress = "pct_mito", verbose = FALSE)
  Seurat_so <- RunPCA(Seurat_so, verbose = FALSE)
  Seurat_so <- RunUMAP(Seurat_so, dims = 1:30, verbose = FALSE)
  
  return(Seurat_so)
}

CellTypeGrapFun <- function(Se_so, iTitle) {
  DimPlot(Se_so, reduction = "umap", group.by = "CellType", shuffle = TRUE) + theme(aspect.ratio = 1) + scale_color_igv() +
  ggtitle(iTitle)
}

PredicRefinedGrapFun <- function(Se_so, iTitle) {
  DimPlot(Se_so, reduction = "umap", group.by = "PredictionRefined", shuffle = TRUE) + theme(aspect.ratio = 1) + 
    scale_color_igv() + ggtitle(iTitle)
}


AML419_dem_tbl <- read_tsv("../Data/hwData/GSM3587950_AML419A-D0.dem.txt.gz", show_col_types = FALSE)
AML419_annot_tbl <- read_tsv("../Data/hwData/GSM3587951_AML419A-D0.anno.txt.gz", show_col_types = FALSE)
Seurat_so = SeuratFun(AML419_dem_tbl, AML419_annot_tbl)
CellTypeGrapFun(Seurat_so, "AML419 sick1 CellType")
PredicRefinedGrapFun(Seurat_so, "AML419 sick1 Prediction")

AML707_dem_tbl <- read_tsv("../Data/hwData/GSM3587969_AML707B-D0.dem.txt.gz", show_col_types = FALSE)
AML707_annot_tbl <- read_tsv("../Data/hwData/GSM3587970_AML707B-D0.anno.txt.gz", show_col_types = FALSE)
Seurat_so = SeuratFun(AML707_dem_tbl, AML707_annot_tbl)
CellTypeGrapFun(Seurat_so, "AML707 sick2 CellType")
PredicRefinedGrapFun(Seurat_so, "AML707 sick2 Prediction")

