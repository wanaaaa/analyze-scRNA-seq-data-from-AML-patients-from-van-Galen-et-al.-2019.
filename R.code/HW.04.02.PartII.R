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

#' ######################################################################
#' ######################################################################
#' Part II
#' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# ======================================================
# Perform sctransform-based analysis
# ======================================================

# Run sctransform.
Seurat_so <- SCTransform(Seurat_so, vars.to.regress = "pct_mito", verbose = FALSE)

# The sctransform results are stored in the new SCT assay.
Seurat_so@assays
DefaultAssay(Seurat_so)

# Perform dimensionality reduction based on the new sctransform normalization.
Seurat_so <- RunPCA(Seurat_so , verbose = FALSE)
Seurat_so <- RunUMAP(Seurat_so , dims = 1:30, verbose = FALSE)

# ======================================================
# Generate sctransform-based UMAPs
# ======================================================
# Check some known population marker genes.
gradient_colors <- c("gray90", "red4")
FeaturePlot(Seurat_so, reduction = "umap", features = "CD34", cols = gradient_colors) +
  theme(aspect.ratio = 1)
FeaturePlot(Seurat_so, reduction = "umap", features = "MPO", cols = gradient_colors) +
  theme(aspect.ratio = 1)
# Check the metadata overlaid onto the UMAP.
FeaturePlot(Seurat_so, reduction = "umap", features = "Score_HSC", cols = gradient_colors) +
  theme(aspect.ratio = 1)
DimPlot(Seurat_so, reduction = "umap", group.by = "CellType", shuffle = TRUE) +
  theme(aspect.ratio = 1) +
  scale_color_igv()
DimPlot(Seurat_so, reduction = "umap", group.by = "PredictionRefined", shuffle = TRUE) +
  theme(aspect.ratio = 1) + scale_color_igv()

# ======================================================
# Get cluster markers (differentially expressed genes)
# ======================================================

# Set the assay to prevent any potential confusion later.
DefaultAssay(Seurat_so) <- "RNA"

# Set the cell identities to a specific grouping
Idents(Seurat_so) <- "CellType"
table(Idents(Seurat_so))

# Find markers for every cluster compared to all remaining cells. There are many parameters you can adjust to make the process more or less stringent.
Seurat_markers <- FindAllMarkers(Seurat_so, only.pos = TRUE, min.pct = 0.2, test.use = "wilcox", logfc.threshold = 0.2, min.cells.group = 10, verbose = FALSE)

# Check the top markers for each cluster (or cell type in this example). Some of the cell types will be missing in the output since they have too few cells.
Seurat_markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)



