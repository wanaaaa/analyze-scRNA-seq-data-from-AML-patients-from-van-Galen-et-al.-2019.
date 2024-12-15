library(readr)

library(tidyverse)
library(Seurat)
library(ggsci)
library(readr)
library(dplyr)

BM3_dem_tbl <- read_tsv("../Data/hwData/GSM3587998_BM3.dem.txt.gz", show_col_types = FALSE)
BM3_dem_tbl[1:8, 1:8]

BM3_annot_tbl <- read_tsv("../Data/hwData/GSM3587999_BM3.anno.txt.gz", show_col_types = FALSE)
head(BM3_annot_tbl)

# Clean up Data
class(BM3_annot_tbl)
BM3_dem_mat <- BM3_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
BM3_dem_mat[1:5, 1:3]
# Check the number of UMIs/counts per cell/column.
BM3_dem_mat %>% colSums() %>% summary()

# The annotation data should be a data frame.
BM3_anno_df <- BM3_annot_tbl %>% column_to_rownames("Cell") %>% as.data.frame()
head(BM3_anno_df)

# Check the values of the annotation columns.
summary(BM3_anno_df$TranscriptomeUMIs)
summary(BM3_anno_df$NumberOfGenes)
table(BM3_anno_df$CellType, useNA = "ifany")
table(BM3_anno_df$PredictionRefined, useNA = "ifany")

# ======================================================
# Create the Seurat object
# ======================================================

# Initialize the Seurat object from the raw counts matrix.
BM3_so <- CreateSeuratObject(counts = BM3_dem_mat, min.cells = 1, min.features = 1)
BM3_so
# By default, all the data is stored in the RNA assay.
BM3_so@assays

# The [[ operator can add columns to object metadata. We can add the percentage of reads that map to the mitochondrial genome.
BM3_so[["pct_mito"]] <- PercentageFeatureSet(BM3_so, pattern = "^MT-")

# The metadata is stored as a data frame in the meta.data slot.
head(BM3_so@meta.data)

# add additional metadata using AddMetaData().
BM3_so <- AddMetaData(BM3_so, metadata = BM3_anno_df)
head(BM3_so@meta.data)

# These cells are already pre-filtered, but we can visualize the QC metrics anyway.
VlnPlot(BM3_so, features = c("nFeature_RNA", "nCount_RNA", "pct_mito"), ncol = 3)

# ======================================================
# Perform standard analysis (log-transformed)
# ======================================================

# The default global-scaling normalization method LogNormalize normalizes the feature expression measurements for each cell
# by the total expression, multiplies by a scale factor of 10,000, and log-transforms the result.
BM3_so <- NormalizeData(BM3_so , normalization.method = "LogNormalize")

# The normalized values are stored in the data slot in the RNA assay
BM3_so@assays$RNA$data[1:8, 1:3]
# Each cell now has 10,000 normalized counts.
BM3_so@assays$RNA$data %>% expm1() %>% colSums() %>% summary()

# Identify highly variable features, scale, and perform dimensionality reduction.
BM3_so <- FindVariableFeatures(BM3_so, verbose = FALSE)
BM3_so <- ScaleData(BM3_so, vars.to.regress = "pct_mito", verbose = FALSE)
BM3_so <- RunPCA(BM3_so, verbose = FALSE)
# BM3_so <- FindNeighbors(BM3_so)
# BM3_so <- FindClusters(BM3_so, resolution = c(0.8, 1, 1.2))
BM3_so <- RunUMAP(BM3_so, dims = 1:30, verbose = FALSE)
BM3_so

# ======================================================
# Perform standard analysis (log-transformed)
# ======================================================

# Check some known population marker genes.
DimPlot(BM3_so, reduction = "umap", group.by = "CellType", shuffle = TRUE) + theme(aspect.ratio = 1) + scale_color_igv()

DimPlot(BM3_so, reduction = "umap", group.by = "PredictionRefined", shuffle = TRUE) +theme(aspect.ratio = 1) + scale_color_igv()

